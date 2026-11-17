/*
 * Poly2Tri Copyright (c) 2009-2022, Poly2Tri Contributors
 * https://github.com/jhasse/poly2tri
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * * Neither the name of Poly2Tri nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "index_codec.h"

#include <algorithm>
#include <cstring>
#include <new>
#include <queue>

namespace p2t {

namespace {

const unsigned char kMagic[4] = { 'P', '2', 'T', 'C' };
const unsigned char kVersion = 1;

const unsigned char kFlagReordered = 1 << 0;

const unsigned char kTagLiteral = 0;
const unsigned char kTagRle = 1;
const unsigned char kTagPattern = 2;

void AppendU64(std::vector<std::uint8_t>& out, std::uint64_t value)
{
  do {
    std::uint8_t byte = static_cast<std::uint8_t>(value & 0x7Fu);
    value >>= 7;
    if (value != 0) {
      byte |= 0x80u;
    }
    out.push_back(byte);
  } while (value != 0);
}

std::uint64_t ZigZagEncode(std::int64_t value)
{
  const std::uint64_t bits = static_cast<std::uint64_t>(value);
  return (bits << 1) ^ (0u - (bits >> 63));
}

std::int64_t ZigZagDecode(std::uint64_t value)
{
  return static_cast<std::int64_t>((value >> 1) ^ (~(value & 1u) + 1u));
}

class ByteReader {
public:
  ByteReader(const std::uint8_t* data, std::size_t size)
    : pointer_(data), end_(data + size)
  {
  }

  std::uint8_t ReadByte()
  {
    if (pointer_ >= end_) {
      throw IndexCodecError("index codec: truncated stream");
    }
    return *pointer_++;
  }

  std::uint64_t ReadU64()
  {
    std::uint64_t result = 0;
    unsigned shift = 0;
    while (true) {
      const std::uint64_t byte = ReadByte();
      if (shift == 63 && (byte & 0x7Eu) != 0) {
        throw IndexCodecError("index codec: varint overflow");
      }
      result |= (byte & 0x7Fu) << shift;
      if ((byte & 0x80u) == 0) {
        return result;
      }
      shift += 7;
      if (shift >= 70) {
        throw IndexCodecError("index codec: varint too long");
      }
    }
  }

  bool AtEnd() const
  {
    return pointer_ == end_;
  }

private:
  const std::uint8_t* pointer_;
  const std::uint8_t* const end_;
};

// Breadth-first numbering of the triangle adjacency graph. Vertices
// sharing a triangle become neighbours; the smallest-id neighbour is
// visited first so the numbering is fully deterministic. Vertices that
// are never referenced keep their relative order in trailing components.
std::vector<std::uint32_t> BuildBfsPermutation(const std::uint32_t* indices,
                                               std::size_t index_count,
                                               std::uint32_t vertex_count)
{
  std::vector<std::vector<std::uint32_t>> adjacency(vertex_count);
  // Only complete triples carry connectivity; a trailing one- or
  // two-index tail is still numbered, just without graph edges.
  for (std::size_t t = 0; t + 3u <= index_count; t += 3) {
    const std::uint32_t a = indices[t];
    const std::uint32_t b = indices[t + 1];
    const std::uint32_t c = indices[t + 2];
    adjacency[a].push_back(b);
    adjacency[a].push_back(c);
    adjacency[b].push_back(a);
    adjacency[b].push_back(c);
    adjacency[c].push_back(a);
    adjacency[c].push_back(b);
  }

  std::vector<std::uint32_t> permutation;
  permutation.reserve(vertex_count);
  std::vector<char> visited(vertex_count, 0);
  std::queue<std::uint32_t> queue;

  for (std::uint32_t start = 0; start < vertex_count; ++start) {
    if (visited[start] != 0) {
      continue;
    }
    visited[start] = 1;
    queue.push(start);
    while (!queue.empty()) {
      const std::uint32_t current = queue.front();
      queue.pop();
      permutation.push_back(current);

      std::vector<std::uint32_t>& neighbours = adjacency[current];
      std::sort(neighbours.begin(), neighbours.end());
      neighbours.erase(std::unique(neighbours.begin(), neighbours.end()),
                       neighbours.end());
      for (const std::uint32_t next : neighbours) {
        if (visited[next] == 0) {
          visited[next] = 1;
          queue.push(next);
        }
      }
    }
  }
  return permutation;
}

} // namespace

IndexCodecResult EncodeIndices(const std::uint32_t* indices,
                               std::size_t index_count,
                               std::uint32_t vertex_count,
                               const IndexCodecOptions& options)
{
  if (vertex_count == 0u && index_count != 0u) {
    throw IndexCodecError("index codec: indices reference an empty vertex set");
  }
  for (std::size_t i = 0; i < index_count; ++i) {
    if (indices[i] >= vertex_count) {
      throw IndexCodecError("index codec: input index out of vertex range");
    }
  }

  IndexCodecResult result;
  std::vector<std::uint8_t>& out = result.bytes;

  const bool reorder = options.reorder_vertices && index_count > 0u;
  std::vector<std::uint32_t> new_id_of_old;
  if (reorder) {
    result.vertex_permutation =
      BuildBfsPermutation(indices, index_count, vertex_count);
    result.reordered = true;

    new_id_of_old.resize(vertex_count);
    for (std::uint32_t new_id = 0; new_id < vertex_count; ++new_id) {
      new_id_of_old[result.vertex_permutation[new_id]] = new_id;
    }
  }

  // Stage 1+2: remap into BFS space and delta-code.
  std::vector<std::int64_t> deltas;
  deltas.reserve(index_count);
  std::int64_t previous = 0;
  for (std::size_t i = 0; i < index_count; ++i) {
    const std::int64_t id =
      reorder ? static_cast<std::int64_t>(new_id_of_old[indices[i]])
              : static_cast<std::int64_t>(indices[i]);
    deltas.push_back(id - previous);
    previous = id;
  }

  // Header.
  for (const unsigned char magic_byte : kMagic) {
    out.push_back(magic_byte);
  }
  out.push_back(kVersion);
  out.push_back(reorder ? kFlagReordered : 0u);
  AppendU64(out, index_count);
  AppendU64(out, vertex_count);
  if (reorder) {
    AppendU64(out, vertex_count);
    for (const std::uint32_t old_id : result.vertex_permutation) {
      AppendU64(out, old_id);
    }
  }

  // Stage 3: greedy run-length encoding over the delta stream.
  // A run is the smallest repeating pattern (length p) that occurs at
  // least twice and covers at least run_threshold entries overall.
  const bool runs_enabled = options.run_threshold >= 2u && index_count > 0u;
  const std::uint64_t threshold =
    runs_enabled ? static_cast<std::uint64_t>(options.run_threshold) : 0u;
  // Pattern search is capped so encode time stays linear in the stream
  // length for large threshold values: a pattern longer than the run
  // threshold could never qualify on its own anyway.
  const std::uint64_t search_cap = std::min<std::uint64_t>(threshold, 64u);

  std::vector<std::int64_t> pending_literals;
  auto flush_literals = [&]() {
    if (pending_literals.empty()) {
      return;
    }
    out.push_back(kTagLiteral);
    out.push_back(static_cast<std::uint8_t>(pending_literals.size()));
    for (const std::int64_t delta : pending_literals) {
      AppendU64(out, ZigZagEncode(delta));
    }
    pending_literals.clear();
  };

  std::size_t position = 0;
  while (position < index_count) {
    std::uint64_t pattern_length = 0;
    std::uint64_t repetitions = 0;

    if (runs_enabled) {
      const std::size_t max_pattern = static_cast<std::size_t>(
        std::min<std::uint64_t>(search_cap, (index_count - position) / 2u));
      for (std::size_t p = 1; p <= max_pattern; ++p) {
        std::size_t r = 1;
        while (position + (r + 1) * p <= index_count) {
          bool matches = true;
          for (std::size_t k = 0; k < p; ++k) {
            if (deltas[position + r * p + k] != deltas[position + k]) {
              matches = false;
              break;
            }
          }
          if (!matches) {
            break;
          }
          ++r;
        }
        if (r >= 2 && static_cast<std::uint64_t>(r) * p >= threshold) {
          pattern_length = p;
          repetitions = r;
          break;
        }
      }
    }

    if (pattern_length == 0u) {
      pending_literals.push_back(deltas[position]);
      ++position;
      if (pending_literals.size() == 255u) {
        flush_literals();
      }
      continue;
    }

    flush_literals();
    const std::size_t total =
      static_cast<std::size_t>(pattern_length * repetitions);
    if (pattern_length == 1u) {
      out.push_back(kTagRle);
      AppendU64(out, repetitions);
    } else {
      out.push_back(kTagPattern);
      AppendU64(out, pattern_length);
      AppendU64(out, repetitions);
    }
    for (std::size_t k = 0; k < pattern_length; ++k) {
      AppendU64(out, ZigZagEncode(deltas[position + k]));
    }
    position += total;
  }
  flush_literals();

  return result;
}

std::vector<std::uint32_t> DecodeIndices(const std::uint8_t* data,
                                         std::size_t size,
                                         std::uint32_t vertex_count)
{
  if (data == nullptr && size != 0u) {
    throw IndexCodecError("index codec: null input buffer");
  }
  ByteReader reader(data, size);

  for (const unsigned char magic_byte : kMagic) {
    if (reader.ReadByte() != magic_byte) {
      throw IndexCodecError("index codec: bad magic");
    }
  }
  if (reader.ReadByte() != kVersion) {
    throw IndexCodecError("index codec: unsupported version");
  }
  const std::uint8_t flags = reader.ReadByte();
  if ((flags & ~kFlagReordered) != 0u) {
    throw IndexCodecError("index codec: unknown header flags");
  }

  const std::uint64_t index_count = reader.ReadU64();
  const std::uint64_t stored_vertex_count = reader.ReadU64();
  if (stored_vertex_count != vertex_count) {
    throw IndexCodecError("index codec: vertex count does not match stream");
  }
  if (index_count > static_cast<std::uint64_t>(-1) / sizeof(std::uint32_t)) {
    throw IndexCodecError("index codec: stream too large");
  }

  std::vector<std::uint32_t> old_id_of_new;
  if ((flags & kFlagReordered) != 0u) {
    const std::uint64_t permutation_count = reader.ReadU64();
    if (permutation_count != vertex_count) {
      throw IndexCodecError("index codec: malformed vertex permutation");
    }
    old_id_of_new.resize(vertex_count);
    std::vector<char> used(vertex_count, 0);
    for (std::uint32_t new_id = 0; new_id < vertex_count; ++new_id) {
      const std::uint64_t old_id = reader.ReadU64();
      if (old_id >= vertex_count || used[static_cast<std::size_t>(old_id)] != 0) {
        throw IndexCodecError("index codec: malformed vertex permutation");
      }
      used[static_cast<std::size_t>(old_id)] = 1;
      old_id_of_new[new_id] = static_cast<std::uint32_t>(old_id);
    }
  }

  std::vector<std::uint32_t> indices;
  try {
    indices.reserve(static_cast<std::size_t>(index_count));
  } catch (const std::bad_alloc&) {
    throw IndexCodecError("index codec: stream too large");
  }

  // Modular arithmetic keeps the running sum well defined even when a
  // corrupt delta would wrap; the range check afterwards rejects it.
  std::uint64_t running = 0;
  std::uint64_t produced = 0;
  auto emit_delta = [&](std::int64_t delta) {
    running = running + static_cast<std::uint64_t>(delta);
    if (running >= vertex_count) {
      throw IndexCodecError("index codec: decoded index out of range");
    }
    indices.push_back(static_cast<std::uint32_t>(running));
    ++produced;
  };

  while (produced < index_count) {
    const std::uint8_t tag = reader.ReadByte();
    if (tag == kTagLiteral) {
      const std::uint64_t count = reader.ReadByte();
      if (count == 0u || produced + count > index_count) {
        throw IndexCodecError("index codec: malformed literal run");
      }
      for (std::uint64_t i = 0; i < count; ++i) {
        emit_delta(ZigZagDecode(reader.ReadU64()));
      }
    } else if (tag == kTagRle || tag == kTagPattern) {
      const std::uint64_t pattern_length =
        (tag == kTagPattern) ? reader.ReadU64() : 1u;
      const std::uint64_t repetitions = reader.ReadU64();
      if (pattern_length == 0u || repetitions < 2u) {
        throw IndexCodecError("index codec: malformed run");
      }
      if (pattern_length > index_count ||
          repetitions > index_count / pattern_length) {
        throw IndexCodecError("index codec: run length overflow");
      }
      const std::uint64_t total = pattern_length * repetitions;
      if (total > index_count - produced) {
        throw IndexCodecError("index codec: run exceeds index count");
      }
      std::int64_t pattern[256];
      if (pattern_length > 256u) {
        throw IndexCodecError("index codec: pattern too long");
      }
      for (std::uint64_t k = 0; k < pattern_length; ++k) {
        pattern[k] = ZigZagDecode(reader.ReadU64());
      }
      for (std::uint64_t i = 0; i < total; ++i) {
        emit_delta(pattern[i % pattern_length]);
      }
    } else {
      throw IndexCodecError("index codec: unknown record tag");
    }
  }

  if (!reader.AtEnd()) {
    throw IndexCodecError("index codec: trailing bytes after stream");
  }

  if ((flags & kFlagReordered) != 0u) {
    for (std::uint32_t& index : indices) {
      index = old_id_of_new[index];
    }
  }
  return indices;
}

} // namespace p2t
