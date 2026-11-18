// Minimal fallback implementing the subset of Boost.Test macros used by the
// Taubin tests, for environments without Boost. When Boost is available the
// test file includes <boost/test/unit_test.hpp> instead and this is unused.
#ifndef MESH_BOOST_TEST_SHIM_H
#define MESH_BOOST_TEST_SHIM_H

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace shim {

struct TestCase {
    std::string name;
    std::function<void()> fn;
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> cases;
    return cases;
}

inline int& failures() {
    static int count = 0;
    return count;
}

struct Registrar {
    Registrar(const std::string& name, std::function<void()> fn) {
        registry().push_back({name, std::move(fn)});
    }
};

struct RequireFailure : std::exception {};

inline void reportFailure(const char* expr, const char* file, int line,
                          bool fatal) {
    ++failures();
    std::cerr << file << ":" << line << ": "
              << (fatal ? "fatal " : "") << "check failed: " << expr << "\n";
    if (fatal) throw RequireFailure();
}

} // namespace shim

#define SHIM_CONCAT_(a, b) a##b
#define SHIM_CONCAT(a, b) SHIM_CONCAT_(a, b)

#define BOOST_AUTO_TEST_CASE(name)                                            \
    static void name();                                                       \
    static ::shim::Registrar SHIM_CONCAT(registrar_, name)(#name, name);     \
    static void name()

#define BOOST_CHECK(expr)                                                     \
    do {                                                                      \
        if (!(expr)) ::shim::reportFailure(#expr, __FILE__, __LINE__, false); \
    } while (0)

#define BOOST_REQUIRE(expr)                                                   \
    do {                                                                      \
        if (!(expr)) ::shim::reportFailure(#expr, __FILE__, __LINE__, true);  \
    } while (0)

#define BOOST_CHECK_SMALL(value, tolerance)                                   \
    do {                                                                      \
        if (std::fabs(value) >= (tolerance))                                  \
            ::shim::reportFailure(#value " ~= 0", __FILE__, __LINE__, false); \
    } while (0)

#define BOOST_REQUIRE_EQUAL(a, b) BOOST_REQUIRE((a) == (b))
#define BOOST_REQUIRE_GE(a, b) BOOST_REQUIRE((a) >= (b))
#define BOOST_REQUIRE_GT(a, b) BOOST_REQUIRE((a) > (b))
#define BOOST_CHECK_EQUAL(a, b) BOOST_CHECK((a) == (b))
#define BOOST_CHECK_GE(a, b) BOOST_CHECK((a) >= (b))
#define BOOST_CHECK_GT(a, b) BOOST_CHECK((a) > (b))
#define BOOST_CHECK_LT(a, b) BOOST_CHECK((a) < (b))

int main() {
    for (const auto& tc : ::shim::registry()) {
        std::cout << "Running " << tc.name << "...\n";
        try {
            tc.fn();
        } catch (const ::shim::RequireFailure&) {
            // failure already recorded
        } catch (const std::exception& e) {
            ++::shim::failures();
            std::cerr << "  uncaught exception: " << e.what() << "\n";
        }
    }
    int failed = ::shim::failures();
    std::cout << ::shim::registry().size() << " test case(s), "
              << failed << " failure(s)\n";
    return failed == 0 ? 0 : 1;
}

#endif
