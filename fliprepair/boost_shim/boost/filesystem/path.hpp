/*
 * Minimal hand-written stand-in for <boost/filesystem/path.hpp>.
 *
 * Implements only the path operations used by unittest/main.cpp:
 * construction from a string, remove_filename(), parent_path(), operator/
 * and string(). Original code, not third-party source.
 */
#pragma once

#include <string>

namespace boost {
namespace filesystem {

class path {
public:
  path() = default;
  path(const char* p) : path_(p) {}
  path(const std::string& p) : path_(p) {}

  path remove_filename() const
  {
    const std::string::size_type pos = path_.find_last_of('/');
    if (pos == std::string::npos) {
      return path();
    }
    return path(path_.substr(0, pos));
  }

  path parent_path() const
  {
    const std::string::size_type pos = path_.find_last_of('/');
    if (pos == std::string::npos) {
      return path();
    }
    if (pos == 0) {
      return path("/");
    }
    return path(path_.substr(0, pos));
  }

  std::string string() const { return path_; }

  path operator/(const path& other) const
  {
    if (path_.empty()) {
      return other;
    }
    if (path_[path_.size() - 1] == '/') {
      return path(path_ + other.path_);
    }
    return path(path_ + "/" + other.path_);
  }

private:
  std::string path_;
};

} // namespace filesystem
} // namespace boost
