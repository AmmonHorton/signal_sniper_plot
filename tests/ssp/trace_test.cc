#include "render/trace.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <set>
#include <vector>

#include "app/compose.h"

namespace ssp {
namespace {

TEST(ReduceTrace, BinsMatchSamplesInEachColumn) {
    std::vector<double> v(10'007);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::sin(0.01 * i) + 0.001 * (i % 13);
    const double xstart = 5.0, xdelta = 0.25;
    Lod lod{Signal(v, xstart, xdelta)};
    const XView view{100.0, 2000.0, 333, 100, CMode::Real, PhaseUnits::Radians};

    // Incremental reduce in uneven steps must equal the one-shot result.
    TraceBins bins, stepped;
    reduce_trace(lod, xstart, xdelta, view, bins);
    stepped.cols.assign(view.width, Bin{});
    TraceReducer r(lod, xstart, xdelta, view, stepped);
    for (int end : {1, 7, 8, 100, 332, 333}) r.run_to(end);
    for (int c = 0; c < view.width; ++c) {
        ASSERT_EQ(stepped.cols[c].count, bins.cols[c].count);
        ASSERT_EQ(stepped.cols[c].lo, bins.cols[c].lo);
        ASSERT_EQ(stepped.cols[c].last, bins.cols[c].last);
    }
    ASSERT_TRUE(stepped.after.valid);

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
    reduce_trace(lod, 0, 1, XView{-0.5, 2.5, 3, 100, CMode::Log20, PhaseUnits::Radians}, bins);
    EXPECT_NEAR(bins.cols[0].hi, 20 * std::log10(5.0), 1e-9);
    EXPECT_NEAR(bins.cols[1].lo, 20 * std::log10(kLogFloor), 1e-9);
    EXPECT_NEAR(bins.cols[2].lo, -20.0, 1e-5);
}

TEST(PaintTrace, ConstantSignalIsOneRow) {
    std::vector<float> v(1000, 0.0f);
    Lod lod{Signal(v)};
    const XView view{0.0, 999.0, 100, 100, CMode::Real, PhaseUnits::Radians};
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
    const XView view{0.4, 0.6, 50, 100, CMode::Real, PhaseUnits::Radians};
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

TEST(ReduceTrace, DotRowsAreExactlyTheRowsHoldingSamples) {
    std::vector<double> v(50'000);
    uint32_t st = 7;
    for (auto& x : v) {
        st = st * 1664525u + 1013904223u;
        x = double(st >> 20) / 4096.0 * 2.0 - 1.0;  // pseudo-random in [-1, 1)
    }
    Lod lod{Signal(v)};
    const XView view{0.0, double(v.size()), 97, 61, CMode::Real, PhaseUnits::Radians};
    TraceBins bins;
    bins.cols.assign(view.width, Bin{});
    bins.occ_words = 1;
    bins.occ.assign(view.width, 0);
    bins.occ_y = {-0.8, 1.3};
    TraceReducer(lod, 0.0, 1.0, view, bins).run_to(view.width);

    std::vector<uint64_t> want(view.width, 0);
    const double col_dx = (view.x1 - view.x0) / view.width;
    for (std::size_t i = 0; i < v.size(); ++i) {
        const int c = static_cast<int>(std::floor(i / col_dx));
        const int r = row_of(v[i], -0.8, 1.3, view.height);
        if (c < view.width && r >= 0 && r < view.height) want[c] |= uint64_t{1} << r;
    }
    for (int c = 0; c < view.width; ++c) EXPECT_EQ(bins.occ[c], want[c]) << "column " << c;
}

/// Distinct pixel rows of the data area holding `color`.
std::set<int> rows_with(const Framebuffer& fb, const Rect& p, uint32_t color) {
    std::set<int> rows;
    for (int y = p.y; y < p.bottom(); ++y) {
        for (int x = p.x; x < p.right(); ++x) {
            if (fb.at(x, y) == color) rows.insert(y);
        }
    }
    return rows;
}

TEST(TraceContent, QpskPhaseDotsShowFourDistinctLevels) {
    // 100k QPSK symbols zoomed all the way out: hundreds of samples per column. Dots must
    // light only the four phase rows; lines would fill the gaps between them.
    std::vector<std::complex<float>> v(100'000);
    uint32_t st = 1;
    for (auto& x : v) {
        st = st * 1664525u + 1013904223u;
        x = std::polar(1.0f, float(M_PI / 4 + (st >> 30) * M_PI / 2));
    }
    PlotOptions o;
    o.cmode = CMode::Phase;
    Signal dots(v);
    dots.style = Style::Dots;
    Screen s(std::make_unique<TraceContent>(std::vector<Signal>{dots}, o));
    Framebuffer fb(500, 300);
    render_headless(fb, s);
    const uint32_t color = Theme{}.trace_color(0);
    EXPECT_EQ(rows_with(fb, s.layout.plot, color).size(), 4u);

    Screen lines(std::make_unique<TraceContent>(std::vector<Signal>{Signal(v)}, o));
    render_headless(fb, lines);
    EXPECT_GT(rows_with(fb, lines.layout.plot, color).size(), 100u);
}

/// Non-empty density cells of trace t.
std::size_t lit_cells(const ReduceResult& r, std::size_t t, uint64_t* total = nullptr) {
    std::size_t n = 0;
    uint64_t sum = 0;
    for (int k = 0; k < r.view.width * r.view.height; ++k) {
        const uint32_t c = r.density[t][k].load();
        n += c != 0;
        sum += c;
    }
    if (total) *total = sum;
    return n;
}

TEST(TraceContent, IrQpskIsFourClusters) {
    std::vector<std::complex<float>> v(200'000);
    uint32_t st = 3;
    for (auto& x : v) {
        st = st * 1664525u + 1013904223u;
        x = std::polar(1.0f, float(M_PI / 4 + (st >> 30) * M_PI / 2));
    }
    PlotOptions o;
    o.cmode = CMode::IR;
    Screen s(std::make_unique<TraceContent>(std::vector<Signal>{Signal(v)}, o));
    Framebuffer fb(500, 300);
    render_headless(fb, s);
    ASSERT_TRUE(s.result->complete());
    uint64_t total = 0;
    const std::size_t cells = lit_cells(*s.result, 0, &total);
    EXPECT_EQ(total, v.size()) << "every sample lands in the autoscaled plane";
    EXPECT_GE(cells, 4u);
    EXPECT_LE(cells, 16u) << "four constellation points, maybe split across pixel edges";
    EXPECT_NEAR(s.xshown.hi, std::cos(M_PI / 4) * 1.04, 0.03) << "x autoscaled to the real part";
}

TEST(TraceContent, IrPlotsOnlyTheChosenTimeWindow) {
    // First half sits at +1, second half at -1: a window in the second half has one cluster.
    std::vector<std::complex<double>> v(10'000);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = {i < 5000 ? 1.0 : -1.0, 0.0};
    PlotOptions o;
    o.cmode = CMode::IR;
    o.xrange = Range{6000, 8000};
    Screen s(std::make_unique<TraceContent>(std::vector<Signal>{Signal(v)}, o));
    Framebuffer fb(400, 300);
    render_headless(fb, s);
    uint64_t total = 0;
    EXPECT_EQ(lit_cells(*s.result, 0, &total), 1u);
    EXPECT_EQ(total, 2001u);
    // Only the -1 samples were used for autoscale (a flat extent widens by 1 each side).
    EXPECT_EQ(s.xshown.lo, -2.0);
    EXPECT_EQ(s.xshown.hi, 0.0);
}

TEST(TraceContent, IrResultsAreReusedOnlyForTheSameWindow) {
    std::vector<std::complex<float>> v(1000, {1.0f, 1.0f});
    PlotOptions o;
    o.cmode = CMode::IR;
    Screen s(std::make_unique<TraceContent>(std::vector<Signal>{Signal(v)}, o));
    Framebuffer fb(300, 200);
    render_headless(fb, s);
    EXPECT_TRUE(s.content->reusable(*s.result, s.request()));
    s.ir[0] = {0, 10};
    EXPECT_FALSE(s.content->reusable(*s.result, s.request()));
}

TEST(TraceContent, AutoscaleCoversEveryVisibleTrace) {
    // Regression: the old code autoscaled from trace 0 only.
    std::vector<float> small(100, 1.0f), big(100, 50.0f);
    PlotOptions o;
    Screen s(std::make_unique<TraceContent>(std::vector<Signal>{Signal(small, 0, 1), Signal(big, 200, 1)}, o));
    Framebuffer fb(400, 300);
    render_headless(fb, s);
    EXPECT_LE(s.views.top().x.lo, 0.0);
    EXPECT_GE(s.views.top().x.hi, 299.0);
    EXPECT_LE(s.yshown.lo, 1.0);
    EXPECT_GE(s.yshown.hi, 50.0);
}

TEST(TraceContent, RejectsBadInput) {
    std::vector<float> v(10);
    PlotOptions o;
    EXPECT_THROW(TraceContent({}, o), std::invalid_argument);
    EXPECT_THROW(TraceContent({Signal(v, 0.0, 0.0)}, o), std::invalid_argument);
    o.yrange = Range{1, 1};
    EXPECT_THROW(TraceContent({Signal(v)}, o), std::invalid_argument);
}

TEST(TraceContent, AutoModePicksMagForComplex) {
    std::vector<std::complex<float>> c(10);
    std::vector<float> r(10);
    EXPECT_EQ(TraceContent({Signal(r)}, {}).initial_cmode(), CMode::Real);
    EXPECT_EQ(TraceContent({Signal(r), Signal(c)}, {}).initial_cmode(), CMode::Mag);
}

}  // namespace
}  // namespace ssp
