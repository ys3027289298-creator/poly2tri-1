/*
 * Minimal hand-written stand-in for <boost/test/unit_test.hpp>.
 *
 * The system has no Boost installed and no network access, so the repo's
 * existing unit tests cannot be built against real Boost. This shim
 * implements exactly the subset of the Boost.Test API used by
 * unittest/main.cpp and unittest/TriangleTest.cpp (test-case registration,
 * BOOST_CHECK/REQUIRE family, a runner main). It is original code written
 * for this purpose, not third-party source. Test files themselves are
 * unmodified.
 */
#pragma once

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace shim_boost_test {

struct RequireFailure : std::exception {
  const char* what() const noexcept override { return "BOOST_REQUIRE failed"; }
};

struct TestCase {
  std::string name;
  void (*fn)();
};

inline std::vector<TestCase>& Registry()
{
  static std::vector<TestCase> registry;
  return registry;
}

inline bool Register(const char* name, void (*fn)())
{
  Registry().push_back(TestCase{name, fn});
  return true;
}

inline int& FailureCount()
{
  static int count = 0;
  return count;
}

inline void RecordFailure(const char* expr, const char* file, int line,
                          const std::string& detail = std::string())
{
  ++FailureCount();
  std::cerr << file << ":" << line << ": check " << expr << " failed";
  if (!detail.empty()) {
    std::cerr << " [" << detail << "]";
  }
  std::cerr << std::endl;
}

template <typename A, typename B>
inline void CheckEqual(const A& a, const B& b, const char* expr_a, const char* expr_b,
                       const char* file, int line, bool required)
{
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif
  if (!(a == b)) {
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    std::ostringstream oss;
    oss << expr_a << " == " << expr_b << " (" << a << " != " << b << ")";
    RecordFailure(oss.str().c_str(), file, line);
    if (required) {
      throw RequireFailure();
    }
  }
}

inline int RunAll(const char* module_name)
{
  int failed_cases = 0;
  for (const auto& test : Registry()) {
    const int failures_before = FailureCount();
    try {
      test.fn();
    } catch (const RequireFailure&) {
      // Already recorded by the failing BOOST_REQUIRE.
    } catch (const std::exception& e) {
      RecordFailure((std::string("unexpected exception: ") + e.what()).c_str(),
                    test.name.c_str(), 0);
    } catch (...) {
      RecordFailure("unexpected unknown exception", test.name.c_str(), 0);
    }
    const bool passed = FailureCount() == failures_before;
    if (!passed) {
      ++failed_cases;
    }
    std::cout << (passed ? "[PASS] " : "[FAIL] ") << test.name << std::endl;
  }
  std::cout << module_name << ": " << (Registry().size() - failed_cases) << "/"
            << Registry().size() << " test cases passed" << std::endl;
  return failed_cases == 0 ? 0 : 1;
}

} // namespace shim_boost_test

#define BOOST_AUTO_TEST_CASE(name)                                              \
  static void name();                                                           \
  static const bool name##_shim_registered =                                    \
      ::shim_boost_test::Register(#name, &name);                                \
  static void name()

#define BOOST_CHECK(cond)                                                       \
  do {                                                                          \
    if (!(cond)) {                                                              \
      ::shim_boost_test::RecordFailure(#cond, __FILE__, __LINE__);              \
    }                                                                           \
  } while (0)

#define BOOST_REQUIRE(cond)                                                     \
  do {                                                                          \
    if (!(cond)) {                                                              \
      ::shim_boost_test::RecordFailure(#cond, __FILE__, __LINE__);              \
      throw ::shim_boost_test::RequireFailure();                                \
    }                                                                           \
  } while (0)

#define BOOST_CHECK_EQUAL(a, b)                                                 \
  ::shim_boost_test::CheckEqual((a), (b), #a, #b, __FILE__, __LINE__, false)

#define BOOST_REQUIRE_EQUAL(a, b)                                               \
  ::shim_boost_test::CheckEqual((a), (b), #a, #b, __FILE__, __LINE__, true)

#define BOOST_CHECK_NO_THROW(expr)                                              \
  do {                                                                          \
    try {                                                                       \
      expr;                                                                     \
    } catch (const std::exception& shim_e) {                                    \
      ::shim_boost_test::RecordFailure(#expr, __FILE__, __LINE__,               \
                                       std::string("threw: ") + shim_e.what()); \
    } catch (...) {                                                             \
      ::shim_boost_test::RecordFailure(#expr, __FILE__, __LINE__,               \
                                       "threw unknown exception");              \
    }                                                                           \
  } while (0)

#define BOOST_CHECK_MESSAGE(cond, msg)                                          \
  do {                                                                          \
    if (!(cond)) {                                                              \
      std::ostringstream shim_oss;                                              \
      shim_oss << msg;                                                          \
      ::shim_boost_test::RecordFailure(#cond, __FILE__, __LINE__,               \
                                       shim_oss.str());                         \
    }                                                                           \
  } while (0)

#ifdef BOOST_TEST_MODULE
#define SHIM_STRINGIZE_IMPL(x) #x
#define SHIM_STRINGIZE(x) SHIM_STRINGIZE_IMPL(x)
int main()
{
  return ::shim_boost_test::RunAll(SHIM_STRINGIZE(BOOST_TEST_MODULE));
}
#endif
