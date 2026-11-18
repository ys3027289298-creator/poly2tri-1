// Minimal boost::filesystem::path-compatible shim.
//
// Used only to run the repository's own unit tests unmodified in an
// environment without Boost (adding third-party dependencies is not
// allowed).  Implements exactly what unittest/main.cpp uses.
#pragma once

#include <string>

namespace boost {
namespace filesystem {

class path {
public:
  path() {}
  path(const char* s) : s_(s) {}
  path(const std::string& s) : s_(s) {}

  path remove_filename() const {
    const std::string::size_type pos = s_.find_last_of('/');
    if (pos == std::string::npos) return path("");
    return path(s_.substr(0, pos + 1));
  }

  path parent_path() const {
    std::string s = s_;
    while (s.size() > 1 && s[s.size() - 1] == '/') s.erase(s.size() - 1);
    const std::string::size_type pos = s.find_last_of('/');
    if (pos == std::string::npos) return path("");
    if (pos == 0) return path("/");
    return path(s.substr(0, pos));
  }

  std::string string() const { return s_; }

  friend path operator/(const path& a, const path& b) {
    if (a.s_.empty() || a.s_[a.s_.size() - 1] == '/') {
      return path(a.s_ + b.s_);
    }
    return path(a.s_ + "/" + b.s_);
  }

private:
  std::string s_;
};

}  // namespace filesystem
}  // namespace boost
