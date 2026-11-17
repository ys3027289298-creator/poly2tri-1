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

#pragma once

#include "dll_symbol.h"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace p2t {

// Topological compression for triangle index buffers.
//
// Pipeline (all stages are lossless):
//   1. Optional BFS vertex reordering so adjacent triangles reference
//      nearby vertex ids, which keeps the subsequent deltas small.
//   2. Delta coding of the index stream (each entry minus its predecessor).
//   3. Run-length encoding: literal runs and repeated delta patterns
//      (covers both constant runs and periodic patterns such as the
//      alternating deltas produced by regular grids).
//   4. LEB128/zig-zag varint packing of every integer.
//
// The byte stream is self-describing and validated on decode: any index
// outside [0, vertex_count), truncated input, trailing garbage, unknown
// flags or a malformed run makes DecodeIndices throw IndexCodecError.

struct IndexCodecOptions {
  // When true, vertices are renumbered along the triangle adjacency graph
  // (breadth-first, smallest-id-first tie breaking) before delta coding.
  // The permutation is stored in the stream and applied backwards by the
  // decoder, so decoded indices are still identical to the input indices.
  bool reorder_vertices = true;

  // Minimum run length, measured in index entries, that is emitted as a
  // run instead of literals. A run needs at least two pattern repetitions
  // and p * repeat_count >= run_threshold entries. Values below 2 disable
  // run-length encoding completely (delta + varint only).
  std::uint32_t run_threshold = 4;
};

struct IndexCodecResult {
  // Compressed byte stream, ready to upload or store.
  std::vector<std::uint8_t> bytes;

  // Empty when IndexCodecOptions::reorder_vertices was false. Otherwise
  // vertex_permutation[new_id] == old_id: reorder the vertex buffer with
  // the same permutation if you want to keep indices in compressed space.
  std::vector<std::uint32_t> vertex_permutation;

  bool reordered = false;
};

class P2T_DLL_SYMBOL IndexCodecError : public std::runtime_error {
public:
  explicit IndexCodecError(const char* message)
    : std::runtime_error(message)
  {
  }
};

// Encodes `index_count` 32-bit indices (typically triangle lists, hence
// a multiple of three, but no topology is assumed) referencing
// `vertex_count` vertices. Throws IndexCodecError on invalid arguments.
P2T_DLL_SYMBOL IndexCodecResult EncodeIndices(
  const std::uint32_t* indices, std::size_t index_count,
  std::uint32_t vertex_count, const IndexCodecOptions& options = IndexCodecOptions{});

// Decodes a stream produced by EncodeIndices and returns the original
// index sequence. `vertex_count` is the bound used for range validation:
// every decoded index is checked against [0, vertex_count). Throws
// IndexCodecError on any malformed input or out-of-range index.
P2T_DLL_SYMBOL std::vector<std::uint32_t> DecodeIndices(
  const std::uint8_t* data, std::size_t size, std::uint32_t vertex_count);

// Convenience overloads for callers holding std::vector data.
inline IndexCodecResult EncodeIndices(const std::vector<std::uint32_t>& indices,
                                      std::uint32_t vertex_count,
                                      const IndexCodecOptions& options = IndexCodecOptions{})
{
  return EncodeIndices(indices.data(), indices.size(), vertex_count, options);
}

inline std::vector<std::uint32_t> DecodeIndices(const std::vector<std::uint8_t>& data,
                                                std::uint32_t vertex_count)
{
  return DecodeIndices(data.data(), data.size(), vertex_count);
}

} // namespace p2t
