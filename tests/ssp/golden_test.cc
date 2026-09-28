// Golden-image tests: render known scenes headlessly and compare pixel-for-pixel.
// Accept new output with:  SSP_UPDATE_GOLDENS=1 bazel run //:ssp_golden_test

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "render/png.h"
#include "app/compose.h"

namespace ssp {
namespace {

constexpr int kW = 480, kH = 300;

void check_golden(const std::string& name, const Framebuffer& fb) {
    const char* ws = std::getenv("BUILD_WORKSPACE_DIRECTORY");
    if (std::getenv("SSP_UPDATE_GOLDENS") && ws) {
        write_png(fb, std::string(ws) + "/tests/golden/" + name + ".png");
        return;
    }
    const std::string path = "tests/golden/" + name + ".png";
    if (!std::ifstream(path)) {
        FAIL() << "missing " << path << "; run SSP_UPDATE_GOLDENS=1 bazel run //:ssp_golden_test";
    }
    const Framebuffer want = read_png(path);
    ASSERT_EQ(want.width(), fb.width());
    ASSERT_EQ(want.height(), fb.height());
    int diffs = 0, fx = -1, fy = -1;
    for (int y = 0; y < fb.height(); ++y) {
        for (int x = 0; x < fb.width(); ++x) {
            if (want.at(x, y) != fb.at(x, y) && diffs++ == 0) { fx = x; fy = y; }
        }
    }
    if (diffs) {
        if (const char* out = std::getenv("TEST_UNDECLARED_OUTPUTS_DIR")) {
            write_png(fb, std::string(out) + "/" + name + ".actual.png");
        }
        FAIL() << name << ": " << diffs << " pixels differ, first at (" << fx << "," << fy << ")";
    }
}

Framebuffer render(std::vector<Signal> sigs, PlotOptions o) {
    Framebuffer fb(kW, kH);
    Screen s(std::make_unique<TraceContent>(std::move(sigs), o), o);
    render_headless(fb, s);
    return fb;
}

TEST(Golden, RealSine) {
    std::vector<float> v(1000);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::sin(2 * M_PI * i / 250.0);
    PlotOptions o;
    o.title = "Real sine";
    check_golden("real_sine", render({Signal(v, 0.0, 1e-3, "sin")}, o));
}

TEST(Golden, TwoComplexTracesDifferentAxes) {
    std::vector<std::complex<float>> a(5000), b(3000);
    for (std::size_t i = 0; i < a.size(); ++i) a[i] = std::polar(1.0f + 0.5f * std::sin(i * 0.004f), i * 0.1f);
    for (std::size_t i = 0; i < b.size(); ++i) b[i] = {0.0f, 2.0f + std::cos(i * 0.01f)};
    PlotOptions o;
    o.title = "Magnitude, xstart/xdelta differ";
    Signal sb(b, 110.0, 0.01, "Imag cos + 2");
    sb.style = Style::Dots;
    check_golden("two_complex_mag", render({Signal(a, 100.0, 0.02, "AM"), sb}, o));
}

TEST(Golden, ZoomedToIndividualSamples) {
    std::vector<double> v(100000);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::sin(i * 0.7) * (i % 5);
    PlotOptions o;
    o.title = "Zoomed in: lines+dots";
    o.xrange = Range{500.5, 540.5};
    Signal s(v, 0.0, 1.0, "v");
    s.style = Style::LinesDots;
    check_golden("zoom_samples", render({s}, o));
}

TEST(Golden, FixedYRangeClipsLines) {
    std::vector<float> v(400);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = 3.0f * std::sin(i * 0.05f);
    PlotOptions o;
    o.title = "yrange -1..1 clips";
    o.yrange = Range{-1, 1};
    o.thickness = 2;
    check_golden("clip_yrange", render({Signal(v)}, o));
}

TEST(Golden, Log20OfToneInNoise) {
    std::vector<std::complex<double>> v(20000);
    uint32_t state = 12345;
    for (std::size_t i = 0; i < v.size(); ++i) {
        state = state * 1664525u + 1013904223u;  // deterministic LCG
        const double noise = (state >> 8) / double(1u << 24) - 0.5;
        v[i] = std::polar(1e-3 * (1.0 + noise) + (i > 8000 && i < 12000 ? 1.0 : 0.0), 0.3 * i);
    }
    PlotOptions o;
    o.title = "20log10 burst";
    o.cmode = CMode::Log20;
    check_golden("log20_burst", render({Signal(v, 0.0, 1e-6)}, o));
}

TEST(Golden, TwoMillionSampleChirp) {
    std::vector<int16_t> v(2'000'000);
    for (std::size_t i = 0; i < v.size(); ++i) {
        const double t = double(i) / v.size();
        v[i] = static_cast<int16_t>(20000 * std::sin(2 * M_PI * (5 + 400 * t) * t) * (0.2 + t));
    }
    PlotOptions o;
    o.title = "2M-sample int16 chirp";
    check_golden("chirp_2m", render({Signal(v, 1e9, 1.0, "chirp")}, o));
}

}  // namespace
}  // namespace ssp
