#include "render/trace.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "render/trace_plot.h"

namespace ssp {
namespace {

TEST(ReduceTrace, BinsMatchSamplesInEachColumn) {
    std::vector<double> v(10'007);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::sin(0.01 * i) + 0.001 * (i % 13);
    const double xstart = 5.0, xdelta = 0.25;
    Lod lod{Signal(v, xstart, xdelta)};
    const XView view{100.0, 2000.0, 333, CMode::Real, PhaseUnits::Radians};

    TraceBins bins;
    std::atomic<int> done{0};
    reduce_trace(lod, xstart, xdelta, view, bins, {}, &done);
    ASSERT_EQ(done.load(), view.width);

    const double col_dx = (view.x1 - view.x0) / view.width;
    std::size_t total = 0;
    for (int c = 0; c < view.width; ++c) {
        Span expect;
        std::size_t count = 0, first = 0, last = 0;
        for (std::size_t i = 0; i < v.size(); ++i) {
            const double x = xstart + i * xdelta;
            if (x >= view.x0 + c * col_dx && x < view.x0 + (c + 1) * col_dx) {
                if (count++ == 0) first = i;
                last = i;
                expect.add(v[i]);
            }
        }
        const Bin& b = bins.cols[c];
        ASSERT_EQ(b.count, count) << "column " << c;
        if (count) {
            EXPECT_EQ(b.lo, expect.lo);
            EXPECT_EQ(b.hi, expect.hi);
            EXPECT_EQ(b.first, v[first]);
            EXPECT_EQ(b.last, v[last]);
        }
        total += count;
    }
    EXPECT_EQ(total, static_cast<std::size_t>((2000.0 - 100.0) / xdelta));
    ASSERT_TRUE(bins.before.valid && bins.after.valid);
    EXPECT_LT(bins.before.x, view.x0);
    EXPECT_GE(bins.after.x, view.x1);
}

TEST(ReduceTrace, LogModesUseMagnitude) {
    std::vector<std::complex<float>> v{{3, 4}, {0, 0}, {0.1f, 0}};
    Lod lod{Signal(v)};
    TraceBins bins;
    reduce_trace(lod, 0, 1, XView{-0.5, 2.5, 3, CMode::Log20, PhaseUnits::Radians}, bins);
    EXPECT_NEAR(bins.cols[0].hi, 20 * std::log10(5.0), 1e-9);
    EXPECT_NEAR(bins.cols[1].lo, 20 * std::log10(kLogFloor), 1e-9);
    EXPECT_NEAR(bins.cols[2].lo, -20.0, 1e-5);
}

TEST(PaintTrace, ConstantSignalIsOneRow) {
    std::vector<float> v(1000, 0.0f);
    Lod lod{Signal(v)};
    const XView view{0.0, 999.0, 100, CMode::Real, PhaseUnits::Radians};
    TraceBins bins;
    reduce_trace(lod, 0, 1, view, bins);
    Framebuffer fb(100, 50);
    paint_trace(fb, bins, 100, view, fb.bounds(), -1.0, 1.0, {Style::Lines, 0xFF0000, 1});
    for (int x = 0; x < 100; ++x) {
        EXPECT_EQ(fb.at(x, 25), 0xFF0000u) << x;
        EXPECT_EQ(fb.at(x, 24), 0u) << x;
        EXPECT_EQ(fb.at(x, 26), 0u) << x;
    }
}

TEST(PaintTrace, ViewBetweenTwoSamplesStillDrawsTheLine) {
    std::vector<float> v{-1.0f, 1.0f};
    Lod lod{Signal(v)};
    const XView view{0.4, 0.6, 50, CMode::Real, PhaseUnits::Radians};
    TraceBins bins;
    reduce_trace(lod, 0, 1, view, bins);
    Framebuffer fb(50, 50);
    paint_trace(fb, bins, 50, view, fb.bounds(), -1.0, 1.0, {Style::Lines, 0xFFFFFF, 1});
    int lit = 0;
    for (int x = 0; x < 50; ++x) {
        for (int y = 0; y < 50; ++y) lit += fb.at(x, y) != 0;
    }
    EXPECT_GE(lit, 50);
}

TEST(TracePlot, AutoscaleCoversEveryVisibleTrace) {
    // Regression: the old code autoscaled from trace 0 only.
    std::vector<float> small(100, 1.0f), big(100, 50.0f);
    TracePlot p({Signal(small, 0, 1), Signal(big, 200, 1)}, {});
    Framebuffer fb(400, 300);
    p.render(fb);
    EXPECT_LE(p.xview().lo, 0.0);
    EXPECT_GE(p.xview().hi, 299.0);
    EXPECT_LE(p.yview().lo, 1.0);
    EXPECT_GE(p.yview().hi, 50.0);
}

TEST(TracePlot, RejectsBadInput) {
    std::vector<float> v(10);
    EXPECT_THROW(TracePlot({}, {}), std::invalid_argument);
    EXPECT_THROW(TracePlot({Signal(v, 0.0, 0.0)}, {}), std::invalid_argument);
    PlotOptions o;
    o.yrange = Range{1, 1};
    EXPECT_THROW(TracePlot({Signal(v)}, o), std::invalid_argument);
}

TEST(TracePlot, AutoModePicksMagForComplex) {
    std::vector<std::complex<float>> c(10);
    std::vector<float> r(10);
    EXPECT_EQ(TracePlot({Signal(r)}, {}).cmode(), CMode::Real);
    EXPECT_EQ(TracePlot({Signal(r), Signal(c)}, {}).cmode(), CMode::Mag);
}

}  // namespace
}  // namespace ssp
