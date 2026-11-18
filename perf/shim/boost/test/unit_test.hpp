// Minimal Boost.Test-compatible shim.
//
// The poly2tri unit tests (unittest/main.cpp, unittest/TriangleTest.cpp)
// are written against Boost.Test, but Boost is not installed in this
// environment and adding third-party dependencies is not allowed.  This
// header implements exactly the macros the test sources use, so the
// repository's own test cases compile and run UNMODIFIED.  It is not a
// third-party library; it is part of this repo's test tooling.
#pragma once

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace boost_shim {

struct TestCase {
  std::string name;
  void (*fn)();
};

inline std::vector<TestCase>& Registry() {
  static std::vector<TestCase> r;
  return r;
}

inline int& CheckFailures() {
  static int f = 0;
  return f;
}

struct Registrar {
  Registrar(const char* name, void (*fn)()) {
    Registry().push_back(TestCase{name, fn});
  }
};

struct RequireFailed {};

inline void Fail(const char* file, int line, const std::string& msg) {
  ++CheckFailures();
  std::cerr << file << "(" << line << "): error: " << msg << std::endl;
}

inline void RequireFail(const char* file, int line, const std::string& msg) {
  Fail(file, line, msg);
  throw RequireFailed();
}

inline int RunAll() {
  int failed_cases = 0;
  std::vector<TestCase>& tests = Registry();
  for (size_t i = 0; i < tests.size(); ++i) {
    const int before = CheckFailures();
    try {
      tests[i].fn();
    } catch (const RequireFailed&) {
    } catch (const std::exception& e) {
      Fail("<uncaught exception>", 0,
           std::string("uncaught std::exception: ") + e.what());
    } catch (...) {
      Fail("<uncaught exception>", 0, "uncaught unknown exception");
    }
    const bool ok = CheckFailures() == before;
    std::cout << (ok ? "[PASS] " : "[FAIL] ") << tests[i].name << std::endl;
    if (!ok) ++failed_cases;
  }
  std::cout << tests.size() << " test case(s), " << failed_cases
            << " failed" << std::endl;
  return failed_cases == 0 ? 0 : 1;
}

}  // namespace boost_shim

#define BOOST_AUTO_TEST_CASE(name)                                     \
  static void name();                                                  \
  static ::boost_shim::Registrar boost_shim_registrar_##name(#name,    \
                                                             &name);   \
  static void name()

#define BOOST_CHECK(expr)                                              \
  do {                                                                 \
    if (!(expr))                                                       \
      ::boost_shim::Fail(__FILE__, __LINE__,                           \
                         "check " #expr " failed");                    \
  } while (0)

#define BOOST_REQUIRE(expr)                                            \
  do {                                                                 \
    if (!(expr))                                                       \
      ::boost_shim::RequireFail(__FILE__, __LINE__,                    \
                                "requirement " #expr " failed");       \
  } while (0)

#define BOOST_CHECK_EQUAL(a, b)                                        \
  do {                                                                 \
    if (!((a) == (b)))                                                 \
      ::boost_shim::Fail(__FILE__, __LINE__,                           \
                         "check " #a " == " #b " failed");             \
  } while (0)

#define BOOST_REQUIRE_EQUAL(a, b)                                      \
  do {                                                                 \
    if (!((a) == (b)))                                                 \
      ::boost_shim::RequireFail(__FILE__, __LINE__,                    \
                                "requirement " #a " == " #b            \
                                " failed");                            \
  } while (0)

#define BOOST_CHECK_NO_THROW(expr)                                     \
  do {                                                                 \
    try {                                                              \
      expr;                                                            \
    } catch (const std::exception& boost_shim_e) {                     \
      ::boost_shim::Fail(__FILE__, __LINE__,                           \
                         std::string("unexpected exception: ") +       \
                             boost_shim_e.what());                     \
    } catch (...) {                                                    \
      ::boost_shim::Fail(__FILE__, __LINE__,                           \
                         "unexpected unknown exception");              \
    }                                                                  \
  } while (0)

#define BOOST_CHECK_MESSAGE(expr, msg)                                 \
  do {                                                                 \
    if (!(expr)) {                                                     \
      std::ostringstream boost_shim_os;                                \
      boost_shim_os << msg;                                            \
      ::boost_shim::Fail(__FILE__, __LINE__,                           \
                         std::string("check " #expr " failed: ") +     \
                             boost_shim_os.str());                     \
    }                                                                  \
  } while (0)

#ifdef BOOST_TEST_MODULE
int main() { return ::boost_shim::RunAll(); }
#endif
