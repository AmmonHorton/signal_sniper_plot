#include "app/raster_content.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "app/compose.h"
#include "app/controller.h"
#include "app/cut.h"

namespace ssp {
namespace {

std::vector<double> noise(std::size_t n, uint32_t seed) {
    std::vector<double> v(n);
    for (auto& x : v) {
        seed = seed * 1664525u + 1013904223u;
        x = double(seed >> 8) / double(1u << 24) * 20.0 - 10.0;
    }
    return v;
}

/// Cells of one axis covered by pixel i, computed independently of the implementation:
/// cells whose centre is inside the pixel, else the cell under the pixel centre.
std::vector<std::size_t> cells_of(int i, double lo, double hi, int n, double start, double delta,
                                  std::size_t count) {
    const double d = (hi - lo) / n, p0 = lo + i * d, p1 = p0 + d;
    std::vector<std::size_t> out;
    for (std::size_t c = 0; c < count; ++c) {
        const double centre = start + (c + 0.5) * delta;
        if (centre >= p0 && centre < p1) out.push_back(c);
    }
    if (out.empty()) {
        const double c = std::floor((lo + (i + 0.5) * d - start) / delta);
        if (c >= 0 && c < double(count)) out.push_back(std::size_t(c));
    }
    return out;
}

/// Reduce `view` of an R x C raster (row-major `v`) and compare every pixel to brute force.
void check_reduce(const std::vector<double>& v, std::size_t cols, Reduce mode, Range x, Range y,
                  int w, int h, CMode cmode = CMode::Real) {
    RasterOptions o;
    o.subsize = cols;
    o.ydelta = 2.0;
    o.reduce = mode;
    o.cmode = cmode;
    Signal sig(v, 5.0, 0.5);
    RasterContent rc(sig, o);
    ReduceRequest q;
    q.view = {x.lo, x.hi, w, h, cmode, PhaseUnits::Radians};
    q.y = y;
    q.reduce = mode;
    auto r = rc.new_result(q);
    rc.reduce(*r, {});
    ASSERT_TRUE(r->complete());

    for (int py = 0; py < h; ++py) {
        const auto rows = cells_of(py, y.lo, y.hi, h, 0.0, 2.0, rc.rows());
        for (int px = 0; px < w; ++px) {
            const auto cs = cells_of(px, x.lo, x.hi, w, 5.0, 0.5, cols);
            double acc = 0;
            int n = 0;
            for (std::size_t rr : rows) {
                for (std::size_t c : cs) {
                    if (rr * cols + c >= v.size()) continue;  // partial last frame
                    double z = display_value(cmode, PhaseUnits::Radians, std::abs(v[rr * cols + c]));
                    if (cmode == CMode::Real) z = v[rr * cols + c];
                    switch (mode) {
                        case Reduce::Max: acc = n ? std::max(acc, z) : z; break;
                        case Reduce::Min: acc = n ? std::min(acc, z) : z; break;
                        case Reduce::MaxAbs: acc = n ? std::max(acc, std::abs(z)) : std::abs(z); break;
                        case Reduce::First: if (!n) acc = z; break;
                        case Reduce::Mean: acc += z; break;
                    }
                    ++n;
                }
            }
            const float got = r->zimg[py * w + px];
            if (n == 0) {
                EXPECT_TRUE(std::isnan(got)) << "pixel " << px << "," << py;
            } else {
                const double want = mode == Reduce::Mean ? acc / n : acc;
                EXPECT_NEAR(got, want, 1e-4 * (1 + std::abs(want))) << "pixel " << px << "," << py;
            }
        }
    }
}

TEST(Raster, ReduceMatchesBruteForceZoomedOut) {
    const std::size_t rows = 37, cols = 53;
    const auto v = noise(rows * cols, 1);
    const Range x{5.0, 5.0 + cols * 0.5}, y{0.0, rows * 2.0};
    for (Reduce m : {Reduce::Max, Reduce::Min, Reduce::Mean, Reduce::MaxAbs, Reduce::First}) {
        SCOPED_TRACE(int(m));
        check_reduce(v, cols, m, x, y, 20, 15);
    }
}

TEST(Raster, ReduceMatchesBruteForceZoomedIn) {
    const auto v = noise(30 * 40, 2);
    check_reduce(v, 40, Reduce::Max, {10.0, 17.0}, {5.0, 21.0}, 90, 70);  // many pixels per cell
}

TEST(Raster, PartialLastFrameIsEmpty) {
    const auto v = noise(10 * 16 + 7, 3);  // 11th frame has 7 of 16 samples
    check_reduce(v, 16, Reduce::Mean, {5.0, 13.0}, {0.0, 22.0}, 16, 11);
}

TEST(Raster, LogModesTransformEachSampleForMean) {
    const auto v = noise(20 * 20, 4);
    check_reduce(v, 20, Reduce::Mean, {5.0, 15.0}, {0.0, 40.0}, 7, 6, CMode::Log20);
    check_reduce(v, 20, Reduce::Max, {5.0, 15.0}, {0.0, 40.0}, 7, 6, CMode::Log20);
}

TEST(Raster, StridedStorageMatchesContiguous) {
    // 12 x 9 cells stored with column stride 2 and 25 samples between frames.
    const std::size_t rows = 12, cols = 9;
    const auto v = noise(rows * cols, 5);
    std::vector<double> padded(rows * 25, -99.0);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) padded[r * 25 + 2 * c] = v[r * cols + c];
    RasterOptions o;
    o.subsize = cols;
    Signal dense(v), sparse(padded);
    sparse.n = rows * cols;
    sparse.stride = 2;
    RasterOptions os = o;
    os.frame_stride = 25;
    RasterContent a(dense, o), b(sparse, os);
    for (std::size_t r = 0; r < rows; ++r) {
        EXPECT_EQ(sample_value(b.row_signal(r, 0, cols), Comp::Re, 4), v[r * cols + 4]);
    }
    EXPECT_EQ(sample_value(b.col_signal(3, 2, rows), Comp::Re, 5), v[7 * cols + 3]);
    ReduceRequest q;
    q.view = {0.0, double(cols), 5, 4, CMode::Real, PhaseUnits::Radians};
    q.y = Range{0.0, double(rows)};
    auto ra = a.new_result(q), rb = b.new_result(q);
    a.reduce(*ra, {});
    b.reduce(*rb, {});
    EXPECT_EQ(ra->zimg, rb->zimg);
}

TEST(Raster, CutSignalsAndCellLookup) {
    const std::size_t rows = 6, cols = 10;
    std::vector<std::complex<float>> v(rows * cols);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = {float(i), 0.0f};
    RasterOptions o;
    o.subsize = cols;
    o.ystart = 100;
    o.ydelta = 10;
    RasterContent rc(Signal(v, 1.0, 0.5), o);
    EXPECT_EQ(rc.rows(), rows);
    const Signal row = rc.row_signal(2, 3, 7);
    EXPECT_EQ(row.n, 4u);
    EXPECT_EQ(row.xstart, 1.0 + 3 * 0.5);
    EXPECT_EQ(sample_value(row, Comp::Re, 0), 23.0);
    const Signal col = rc.col_signal(4, 1, 5);
    EXPECT_EQ(col.n, 4u);
    EXPECT_EQ(col.xstart, 110.0);
    EXPECT_EQ(col.xdelta, 10.0);
    EXPECT_EQ(sample_value(col, Comp::Re, 2), 34.0);

    EXPECT_EQ(rc.col_at(1.0 + 4.75 * 0.5, false), std::optional<std::size_t>(4));
    EXPECT_FALSE(rc.col_at(0.9, false));
    EXPECT_EQ(rc.row_at(125.0, false), std::optional<std::size_t>(2));
    EXPECT_EQ(rc.cols_in({1.9, 2.6}, false), (SampleRange{1, 4}));  // cells 1..3 overlap
    EXPECT_EQ(rc.rows_in({0.0, 1e9}, false), (SampleRange{0, rows}));
    Settings s = rc.initial_settings();
    EXPECT_EQ(s.cmode, CMode::Mag);
    EXPECT_NEAR(*rc.z_at(1.0 + 7.2 * 0.5, 131.0, s), 37.0, 1e-9);
}

TEST(Raster, HomeCoversAllCellsFrameZeroAtTop) {
    std::vector<float> v(4 * 8);
    RasterOptions o;
    o.subsize = 8;
    o.ystart = 50;
    o.ydelta = 5;
    Screen s(std::make_unique<RasterContent>(Signal(v, 0.0, 2.0), o));
    EXPECT_EQ(s.views.top().x.hi, 16.0);
    ASSERT_TRUE(s.views.top().y);
    EXPECT_EQ(s.views.top().y->lo, 50.0);
    EXPECT_EQ(s.views.top().y->hi, 70.0);
    EXPECT_TRUE(s.ydown());
}

TEST(Raster, RejectsBadOptions) {
    std::vector<float> v(10);
    RasterOptions o;
    EXPECT_THROW(RasterContent(Signal(v), o), std::invalid_argument);  // no subsize
    o.subsize = 5;
    o.cmode = CMode::IR;
    EXPECT_THROW(RasterContent(Signal(v), o), std::invalid_argument);
    o.cmode = CMode::Auto;
    o.ydelta = 0;
    EXPECT_THROW(RasterContent(Signal(v), o), std::invalid_argument);
}

/// A 6-row x 10-column raster (value = 100*row + col), rendered, then zoomed to columns
/// 2..6 and rows 1..4.
class CutTest : public ::testing::Test {
protected:
    void SetUp() override {
        for (std::size_t i = 0; i < v_.size(); ++i) v_[i] = float(100 * (i / 10) + i % 10);
        RasterOptions o;
        o.title = "R";
        o.subsize = 10;
        o.ystart = 50;
        o.ydelta = 5;
        s_ = std::make_unique<Screen>(std::make_unique<RasterContent>(Signal(v_, 0.0, 1.0), o));
        s_->views.push(View{{2.0, 7.0}, Range{55.0, 75.0}});
        render_headless(fb_, *s_);
    }
    std::vector<float> v_ = std::vector<float>(57);  // last frame has 7 of 10 columns
    Framebuffer fb_{400, 300};
    std::unique_ptr<Screen> s_;
};

TEST_F(CutTest, XCutIsTheRowInsideTheZoomBox) {
    auto cut = make_cut(*s_, true, 3);
    const Signal& sig = cut->traces()->signal(0);
    ASSERT_EQ(sig.n, 5u);  // columns 2..6
    EXPECT_EQ(sample_value(sig, Comp::Re, 0), 302.0);
    EXPECT_EQ(sig.xstart, 2.0);
    EXPECT_NE(cut->set.title.find("x-cut row 3"), std::string::npos);
    EXPECT_EQ(sig.name, "row 3");
    EXPECT_FALSE(cut->set.legend) << "hidden at first, 'l' shows it";
    EXPECT_EQ(cut->views.level(), 0u) << "the cut's home is the box";
    EXPECT_EQ(cut->views.top().x.lo, 2.0);
    EXPECT_EQ(cut->views.top().x.hi, 6.0);
}

TEST_F(CutTest, YCutIsTheColumnInsideTheZoomBoxAndSkipsMissingCells) {
    auto cut = make_cut(*s_, false, 8);  // rows 1..4 in the box; the partial row 5 isn't
    const Signal& sig = cut->traces()->signal(0);
    ASSERT_EQ(sig.n, 4u);
    EXPECT_EQ(sample_value(sig, Comp::Re, 2), 308.0);
    EXPECT_EQ(sig.name, "column 8");
    EXPECT_EQ(sig.xstart, 55.0);
    EXPECT_EQ(sig.xdelta, 5.0);

    s_->views.push(View{{0.0, 10.0}, Range{50.0, 80.0}});  // all six rows
    render_headless(fb_, *s_);
    EXPECT_EQ(make_cut(*s_, false, 8)->traces()->signal(0).n, 5u) << "row 5 has no column 8";
    EXPECT_EQ(make_cut(*s_, false, 3)->traces()->signal(0).n, 6u);
}

TEST_F(CutTest, SteppingStaysInsideTheBox) {
    auto cut = make_cut(*s_, true, 3);
    cut->views.push(View{{3.0, 5.0}, std::nullopt});
    EXPECT_TRUE(step_cut(*cut, 1));
    EXPECT_EQ(sample_value(cut->traces()->signal(0), Comp::Re, 0), 402.0);
    EXPECT_EQ(cut->traces()->signal(0).name, "row 4");
    EXPECT_EQ(cut->views.level(), 1u) << "stepping keeps the cut's zoom";
    EXPECT_FALSE(step_cut(*cut, 1)) << "row 5 is outside the box";
    EXPECT_FALSE(cut->message.empty());
    EXPECT_TRUE(step_cut(*cut, -3));
    EXPECT_FALSE(step_cut(*cut, -1));
}

TEST_F(CutTest, KeysRequestCutsUnderThePointer) {
    Screen& s = *s_;
    const Rect p = s.layout.plot;
    InputEvent m{InputEvent::Type::Motion, p.x + p.w / 2, p.y + 1};
    handle_event(s, m, false);
    InputEvent x{InputEvent::Type::Key};
    x.key = 'x';
    EXPECT_TRUE(handle_event(s, x, false) & kCut);
    ASSERT_TRUE(s.cut_request);
    EXPECT_TRUE(s.cut_request->first);
    EXPECT_EQ(s.cut_request->second, 1u) << "top of the box is row 1 (frame 0 is at the top)";
    s.cut_request.reset();
    InputEvent y = x;
    y.key = 'y';
    handle_event(s, y, false);
    EXPECT_EQ(s.cut_request->second, 4u);  // middle of columns 2..6

    handle_event(s, InputEvent{InputEvent::Type::Leave}, false);
    s.cut_request.reset();
    EXPECT_EQ(handle_event(s, x, false) & kCut, 0u) << "no pointer, no cut";
}

TEST_F(CutTest, ColormapZWindowAndModes) {
    Screen& s = *s_;
    InputEvent k{InputEvent::Type::Key};
    k.key = 'c';
    EXPECT_TRUE(handle_event(s, k, false) & kRepaint);
    EXPECT_EQ(s.set.cmap, Colormap::ColorWheel);
    const Range z = s.zshown;
    k.key = ']';
    handle_event(s, k, false);
    ASSERT_TRUE(s.set.zfixed);
    EXPECT_NEAR(s.set.zfixed->lo, z.lo + 0.1 * (z.hi - z.lo), 1e-9);
    k.key = '5';
    EXPECT_EQ(handle_event(s, k, false) & kReduce, 0u) << "no IR for rasters";
    EXPECT_EQ(s.set.cmode, CMode::Real);
    action::open_prompt(s, PromptState::Kind::ZRange);
    for (char c : std::string("-3 9")) {
        k.key = static_cast<uint8_t>(c);
        handle_event(s, k, false);
    }
    k.key = 0xff0d;
    handle_event(s, k, false);
    EXPECT_EQ(s.set.zfixed->lo, -3.0);
    EXPECT_EQ(s.set.zfixed->hi, 9.0);
}

TEST_F(CutTest, BoxZoomOnARasterUsesYDown) {
    Screen& s = *s_;
    const Rect p = s.layout.plot;
    handle_event(s, InputEvent{InputEvent::Type::Press, p.x + 10, p.y + 5, 1}, false);
    EXPECT_TRUE(handle_event(s, InputEvent{InputEvent::Type::Release, p.x + 100, p.y + p.h / 2, 1}, false) & kReduce);
    const Range y = *s.views.top().y;
    EXPECT_NEAR(y.lo, 55.0 + 20.0 * 5 / p.h, 1e-9) << "top of the screen is the low y end";
    EXPECT_NEAR(y.hi, 55.0 + 20.0 * (p.h / 2) / p.h, 1e-9);
}

}  // namespace
}  // namespace ssp
