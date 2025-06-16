#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <complex>
#include "plot.h"

TEST(Plot, two_sinusoids) {
    constexpr size_t N = 1e4;
    std::vector<std::complex<float>> signal;
    for (size_t i = 0; i < N; ++i) {
        float t = i * 0.01f;
        signal.push_back({std::sin(2 * M_PI * t), 0.0f});
    }

    for (size_t i = 0; i < N; ++i) {
        float t = i * 0.01f;
        signal.push_back({0.0f, std::sin(2 * M_PI * t) + 1});
    }


    xplot::plot_buffer(
        signal.data(),
        signal.size(),
        sizeof(float),
        true,  // not complex
        true,   // is float
        100.0,  // J1950 start time (2020)
        0.01,           // xdelta (s)
        "Signal Plot - Frame 37",
        2,
        std::make_optional<std::pair<double, double>>({-5, 5}),
        2
    );
}