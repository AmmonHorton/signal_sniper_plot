#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <complex>
#include "plot.h"

TEST(Plot, two_sinusoids) {
    constexpr size_t N = 1e4;
    constexpr double xstart1 = 100.0;
    constexpr double xdelta1 = 0.02;

    std::vector<std::complex<float>> signal1, signal2;
    for (size_t i = 0; i < N; ++i) {
        float t = i * xdelta1;
        signal1.emplace_back(std::sin(2 * M_PI * t), 0.0f);
    }

    // Make a signal with a different start and smaller xdelta

    constexpr double xstart2 = 110;
    constexpr double xdelta2 = 0.01;

    for (size_t i = 0; i < N; ++i) {
        float t = i * xdelta2;
        signal2.emplace_back(0.0f, std::sin(2 * M_PI * t) + 1.0f);
    }

    std::vector<xplot::Trace> traces;

    // Convert to Trace using process_buffer
    traces.push_back(xplot::process_buffer(signal1.data(), N, sizeof(std::complex<float>),
                                           true, true, xstart1, xdelta1));
    traces.back().label = "Real Sin";

    traces.push_back(xplot::process_buffer(signal2.data(), N, sizeof(std::complex<float>),
                                           true, true, xstart2, xdelta2));
    traces.back().label = "Imag Sin + 1";
    traces.back().style = xplot::PlotRenderStyle::DOTS;

    xplot::plot_buffer_traces(
        traces,
        "Signal Plot - Frame 37",
        2,
        std::make_optional<std::pair<double, double>>({-5, 5}),  // y-range
        std::nullopt,  // x-range
        xplot::XAxisMode::TIME,
        xplot::PlotMode::Magnitude
    );
}
