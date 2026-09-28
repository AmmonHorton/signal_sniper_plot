// Rough timing of the headless xplot pipeline on a large complex signal.
//   bazel run -c opt //:ssp_bench -- [num_samples]
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "render/trace_plot.h"

using Clock = std::chrono::steady_clock;

static double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

int main(int argc, char** argv) {
    const std::size_t n = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 200'000'000;
    std::vector<std::complex<float>> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = std::polar(1.0f + 0.1f * (i % 97), 0.001f * (i % 6283));
    std::printf("%zu complex64 samples (%.1f GB)\n", n, n * 8 / 1e9);

    ssp::Framebuffer fb(1920, 1080);
    // One TracePlot kept across renders so its pyramid is reused, as the interactive app will.
    ssp::TracePlot p({ssp::Signal(v)}, {});
    auto t0 = Clock::now();
    p.render(fb);
    std::printf("  %-44s %9.1f ms\n", "first full view (scan + build pyramid)", ms_since(t0));
    t0 = Clock::now();
    p.render(fb);
    std::printf("  %-44s %9.1f ms\n", "full view again (pyramid)", ms_since(t0));

    const ssp::XView full{0.0, double(n - 1), 1920 - 88, ssp::CMode::Mag, ssp::PhaseUnits::Radians};
    for (double frac : {0.5, 0.01, 1e-4, 1e-6}) {
        ssp::XView z = full;
        z.x0 = 0.3 * n;
        z.x1 = z.x0 + frac * n;
        t0 = Clock::now();
        p.reduce(z);
        char label[64];
        std::snprintf(label, sizeof label, "reduce zoom to %g of data", frac);
        std::printf("  %-44s %9.1f ms\n", label, ms_since(t0));
    }
}
