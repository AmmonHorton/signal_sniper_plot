// Rough timing of the headless xplot pipeline on a large complex signal.
//   bazel run -c opt //:ssp_bench -- [num_samples]
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "app/compose.h"
#include "app/raster_content.h"

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
    // One Screen kept across renders so its pyramid is reused, as the interactive app does.
    ssp::PlotOptions o;
    ssp::Screen s(std::make_unique<ssp::TraceContent>(std::vector<ssp::Signal>{ssp::Signal(v)}, o));
    auto t0 = Clock::now();
    ssp::render_headless(fb, s);
    std::printf("  %-44s %9.1f ms\n", "first full view (scan + build pyramid)", ms_since(t0));
    t0 = Clock::now();
    ssp::render_headless(fb, s);
    std::printf("  %-44s %9.1f ms\n", "full view again (pyramid)", ms_since(t0));

    ssp::XView z = s.xview();
    for (double frac : {0.5, 0.01, 1e-4, 1e-6}) {
        z.x0 = 0.3 * n;
        z.x1 = z.x0 + frac * n;
        ssp::ReduceRequest rq;
        rq.view = z;
        auto r = s.content->new_result(rq);
        t0 = Clock::now();
        s.content->reduce(*r, {});
        char label[64];
        std::snprintf(label, sizeof label, "reduce zoom to %g of data", frac);
        std::printf("  %-44s %9.1f ms\n", label, ms_since(t0));
    }

    // Dots mode reads every sample in view to light only rows that hold samples.
    uint32_t st = 1;
    for (auto& x : v) {
        st = st * 1664525u + 1013904223u;
        x = std::polar(1.0f, float(M_PI / 4 + (st >> 30) * M_PI / 2));  // QPSK
    }
    {
        ssp::PlotOptions po;
        po.cmode = ssp::CMode::IR;
        ssp::Screen q(std::make_unique<ssp::TraceContent>(std::vector<ssp::Signal>{ssp::Signal(v)}, po));
        for (const char* when : {"IR first view (whole signal)", "IR again"}) {
            t0 = Clock::now();
            ssp::render_headless(fb, q);
            std::printf("  %-44s %9.1f ms\n", when, ms_since(t0));
        }
    }
    {
        // The same samples as a raster of 4096-sample frames (~49k frames).
        ssp::RasterOptions ro;
        ro.subsize = 4096;
        ro.cmode = ssp::CMode::Log20;
        ssp::Screen q(std::make_unique<ssp::RasterContent>(ssp::Signal(v), ro));
        for (const char* when : {"raster 20log, whole (max)", "raster again"}) {
            t0 = Clock::now();
            ssp::render_headless(fb, q);
            std::printf("  %-44s %9.1f ms\n", when, ms_since(t0));
        }
        q.views.push(ssp::View{{0.0, 400.0}, ssp::Range{1000.0, 2000.0}});
        t0 = Clock::now();
        ssp::render_headless(fb, q);
        std::printf("  %-44s %9.1f ms\n", "raster zoomed to 400 x 1000 cells", ms_since(t0));
    }
    std::printf("QPSK, phase mode:\n");
    for (ssp::Style style : {ssp::Style::Lines, ssp::Style::Dots}) {
        ssp::Signal sig(v);
        sig.style = style;
        ssp::PlotOptions po;
        po.cmode = ssp::CMode::Phase;
        ssp::Screen q(std::make_unique<ssp::TraceContent>(std::vector<ssp::Signal>{sig}, po));
        const char* name = style == ssp::Style::Dots ? "dots " : "lines";
        for (const char* when : {"first view", "again"}) {
            t0 = Clock::now();
            ssp::render_headless(fb, q);
            std::printf("  %s %-38s %9.1f ms\n", name, when, ms_since(t0));
        }
        for (double frac : {0.01, 1e-4}) {
            ssp::XView zz = q.xview();
            zz.x0 = 0.3 * n;
            zz.x1 = zz.x0 + frac * n;
            ssp::ReduceRequest rq;
            rq.view = zz;
            auto r = q.content->new_result(rq);
            t0 = Clock::now();
            q.content->reduce(*r, {});
            std::printf("  %s zoom to %-30g %9.1f ms\n", name, frac, ms_since(t0));
        }
    }
}
