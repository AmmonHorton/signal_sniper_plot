#include "ssp/core/lod.h"

#include <gtest/gtest.h>

#include <atomic>
#include <cmath>
#include <complex>
#include <limits>
#include <random>
#include <vector>

namespace ssp {
namespace {

/// Straightforward reference: component of every sample, one at a time.
Span brute(const Signal& s, Comp c, std::size_t i0, std::size_t i1) {
    Span out;
    for (std::size_t i = i0; i < std::min(i1, s.n); ++i) out.add(sample_value(s, c, i));
    return out;
}

void expect_same(const Span& a, const Span& b, const std::string& what) {
    ASSERT_EQ(a.empty(), b.empty()) << what;
    if (a.empty()) return;
    EXPECT_EQ(a.lo, b.lo) << what;
    EXPECT_EQ(a.hi, b.hi) << what;
}

/// Random queries (short, block-straddling, and whole-signal) for every component, both
/// while the pyramid is partial and after it is complete.
void check_signal(const Signal& s, uint32_t seed) {
    Lod lod(s);
    std::mt19937 rng(seed);
    auto rand_below = [&](std::size_t n) { return n == 0 ? 0 : rng() % (n + 1); };
    const Comp comps[] = {Comp::Re, Comp::Im, Comp::Mag, Comp::Phase};
    for (int round = 0; round < 2; ++round) {
        for (Comp c : comps) {
            for (int q = 0; q < 60; ++q) {
                std::size_t a = rand_below(s.n), b = rand_below(s.n);
                if (q % 3 == 0) b = std::min(s.n, a + rng() % 50);  // short range
                if (a > b) std::swap(a, b);
                expect_same(lod.query(c, a, b), brute(s, c, a, b),
                            "comp " + std::to_string(int(c)) + " [" + std::to_string(a) + "," +
                                std::to_string(b) + ") round " + std::to_string(round));
            }
            // Whole signal completes the pyramid (first round) and uses the upper levels (second).
            expect_same(lod.query(c, 0, s.n), brute(s, c, 0, s.n), "whole signal");
        }
    }
    if (s.n > 0) EXPECT_EQ(lod.coverage(s.complex ? Comp::Mag : Comp::Re), 1.0);
}

template <class T>
std::vector<T> ramp_noise(std::size_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<T> v(n);
    for (auto& x : v) {
        const double u = static_cast<double>(rng()) / std::numeric_limits<uint32_t>::max();
        if constexpr (std::is_floating_point_v<T>) {
            x = static_cast<T>(200.0 * (u - 0.5));
        } else {
            x = static_cast<T>(std::numeric_limits<T>::min() +
                               u * (static_cast<double>(std::numeric_limits<T>::max()) -
                                    std::numeric_limits<T>::min()));
        }
    }
    return v;
}

TEST(Lod, MatchesBruteForceForEveryRealDtype) {
    constexpr std::size_t N = 50'000;  // not a multiple of the block size
    auto i8 = ramp_noise<int8_t>(N, 1);
    auto u8 = ramp_noise<uint8_t>(N, 2);
    auto i16 = ramp_noise<int16_t>(N, 3);
    auto u16 = ramp_noise<uint16_t>(N, 4);
    auto i32 = ramp_noise<int32_t>(N, 5);
    auto u32 = ramp_noise<uint32_t>(N, 6);
    auto i64 = ramp_noise<int64_t>(N, 7);
    auto u64 = ramp_noise<uint64_t>(N, 8);
    auto f32 = ramp_noise<float>(N, 9);
    auto f64 = ramp_noise<double>(N, 10);
    check_signal(Signal(i8), 11);
    check_signal(Signal(u8), 12);
    check_signal(Signal(i16), 13);
    check_signal(Signal(u16), 14);
    check_signal(Signal(i32), 15);
    check_signal(Signal(u32), 16);
    check_signal(Signal(i64), 17);
    check_signal(Signal(u64), 18);
    check_signal(Signal(f32), 19);
    check_signal(Signal(f64), 20);
}

TEST(Lod, MatchesBruteForceForComplex) {
    constexpr std::size_t N = 40'000;
    auto f = ramp_noise<float>(2 * N, 21);
    auto s16 = ramp_noise<int16_t>(2 * N, 22);
    check_signal(Signal(reinterpret_cast<const std::complex<float>*>(f.data()), N), 23);
    check_signal(Signal(s16.data(), DType::I16, true, N), 24);
}

TEST(Lod, StridedViewIsAColumnOfAMatrix) {
    constexpr std::size_t rows = 5000, cols = 7;
    auto m = ramp_noise<double>(rows * cols, 30);
    Signal col(m.data() + 3, rows);
    col.stride = cols;
    check_signal(col, 31);
    EXPECT_EQ(sample_value(col, Comp::Re, 10), m[10 * cols + 3]);
}

TEST(Lod, IgnoresNaNAndInf) {
    std::vector<float> v(3000, 1.0f);
    v[5] = std::nanf("");
    v[2000] = std::numeric_limits<float>::infinity();
    v[2500] = -3.0f;
    Lod lod{Signal(v)};
    const Span s = lod.query(Comp::Re, 0, v.size());
    EXPECT_EQ(s.lo, -3.0);
    EXPECT_EQ(s.hi, 1.0);
}

TEST(Lod, SmallAndEmptySignals) {
    std::vector<double> one{4.0};
    Lod lod{Signal(one)};
    EXPECT_EQ(lod.query(Comp::Re, 0, 1).lo, 4.0);
    EXPECT_TRUE(lod.query(Comp::Re, 1, 1).empty());
    Lod empty{Signal(one.data(), 0)};
    EXPECT_TRUE(empty.query(Comp::Mag, 0, 10).empty());
}

TEST(Lod, CancelThrowsAndLeavesPyramidConsistent) {
    auto v = ramp_noise<float>(300'000, 40);
    Lod lod{Signal(v)};
    std::atomic<uint64_t> gen{2};
    const CancelToken stale(gen, 1);  // already superseded
    EXPECT_THROW(lod.query(Comp::Re, 0, v.size(), stale), Cancelled);
    expect_same(lod.query(Comp::Re, 0, v.size()), brute(Signal(v), Comp::Re, 0, v.size()),
                "after cancel");
}

TEST(Lod, PhaseUsesAtan2ImagReal) {
    // Regression: the old code computed atan2(real, imag).
    std::vector<std::complex<float>> v{{0.0f, 1.0f}};
    EXPECT_DOUBLE_EQ(sample_value(Signal(v), Comp::Phase, 0), M_PI / 2);
}

TEST(Component, FastAtan2IsWithinTwoMicroradians) {
    double worst = 0.0;
    for (int i = 0; i <= 100'000; ++i) {
        const double t = -M_PI + 2 * M_PI * i / 100'000;
        for (double r : {1e-6, 1.0, 1e9}) {
            const double x = r * std::cos(t), y = r * std::sin(t);
            worst = std::max(worst, std::abs(fast_atan2(y, x) - std::atan2(y, x)));
        }
    }
    EXPECT_LT(worst, 2e-6);
    EXPECT_EQ(fast_atan2(0.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(fast_atan2(0.0, -1.0), M_PI);
    EXPECT_DOUBLE_EQ(fast_atan2(-1.0, 0.0), -M_PI / 2);
}

}  // namespace
}  // namespace ssp
