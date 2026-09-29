// Drives the real event loop (run_plot) with a fake window system.

#include "ssp/app/session.h"

#include <gtest/gtest.h>
#include <poll.h>

#include <chrono>
#include <complex>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <thread>
#include <vector>

#include "ssp/app/compose.h"
#include "ssp/app/job_runner.h"
#include "ssp/app/raster_content.h"
#include "ssp/render/png.h"

namespace ssp {
namespace {

using T = InputEvent::Type;

/// Hands the loop one scripted batch of events each time it goes idle (waits with no
/// timeout), and records what was on screen at that moment. Closes the window when the
/// script runs out.
class FakeBackend : public Backend {
public:
    FakeBackend(int w, int h) : w_(w), h_(h) {}

    std::vector<std::vector<InputEvent>> script;
    std::vector<InputEvent> while_busy;  ///< Delivered once, the first time the loop is busy.
    std::vector<Framebuffer> idle_frames;
    int busy_waits = 0;

    int width() const override { return w_; }
    int height() const override { return h_; }

    void wait_events(std::vector<InputEvent>& out, int timeout_ms, int extra_fd) override {
        if (timeout_ms < 0) {
            idle_frames.push_back(screen_);
            if (next_ < script.size()) {
                out = script[next_++];
            } else {
                out.push_back({T::Close});
            }
            return;
        }
        ++busy_waits;
        if (!while_busy.empty()) {
            out = std::move(while_busy);
            while_busy.clear();
            return;
        }
        pollfd fd{extra_fd, POLLIN, 0};
        ::poll(&fd, 1, timeout_ms);
    }

    void present(const Framebuffer& fb, const std::vector<Rect>&) override { screen_ = fb; }

private:
    int w_, h_;
    std::size_t next_ = 0;
    Framebuffer screen_;
};

/// Pixels of the data area only (the readout/status rows legitimately differ).
bool same_plot_area(const Framebuffer& a, const Framebuffer& b, const Rect& p) {
    for (int y = p.y; y < p.bottom(); ++y) {
        for (int x = p.x; x < p.right(); ++x) {
            if (a.at(x, y) != b.at(x, y)) return false;
        }
    }
    return true;
}

Framebuffer headless(const std::vector<Signal>& sigs, const PlotOptions& o, int w, int h) {
    Screen s(std::make_unique<TraceContent>(sigs, o));
    Framebuffer fb(w, h);
    render_headless(fb, s);
    return fb;
}

constexpr int kW = 500, kH = 320;

std::vector<double> sine(std::size_t n) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = std::sin(i * 0.01) * (1 + i % 7);
    return v;
}

TEST(Session, InteractiveFrameMatchesHeadlessRender) {
    const auto v = sine(100'000);
    const std::vector<Signal> sigs{Signal(v, 0.0, 1.0, "s")};
    FakeBackend be(kW, kH);
    run_plot(sigs, {}, be);
    ASSERT_EQ(be.idle_frames.size(), 1u);
    const Rect p = compute_layout(kW, kH).plot;
    EXPECT_TRUE(same_plot_area(be.idle_frames[0], headless(sigs, {}, kW, kH), p));
}

TEST(Session, ZoomThenUnzoomRestoresThePicture) {
    const auto v = sine(100'000);
    const std::vector<Signal> sigs{Signal(v)};
    const Rect p = compute_layout(kW, kH).plot;
    FakeBackend be(kW, kH);
    be.script = {
        {{T::Press, p.x + 100, p.y + 20, 1}, {T::Release, p.x + 200, p.y + 150, 1}, {T::Leave}},
        {{T::Press, p.x + 100, p.y + 20, 3}, {T::Release, p.x + 100, p.y + 20, 3}, {T::Leave}},
    };
    run_plot(sigs, {}, be);
    ASSERT_EQ(be.idle_frames.size(), 3u);
    EXPECT_FALSE(same_plot_area(be.idle_frames[0], be.idle_frames[1], p));
    EXPECT_TRUE(same_plot_area(be.idle_frames[0], be.idle_frames[2], p));
}

TEST(Session, UnzoomReusesTheCachedLevelWithoutRecomputing) {
    // Dots mode reads every sample on each reduce, so a recompute always keeps the loop busy.
    std::vector<std::complex<float>> v(3'000'000);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::polar(1.0f, float(i % 4) * 1.5708f);
    Signal sig(v);
    sig.style = Style::Dots;
    PlotOptions o;
    o.cmode = CMode::Phase;
    const Rect p = compute_layout(kW, kH).plot;

    class Counting : public FakeBackend {
    public:
        using FakeBackend::FakeBackend;
        std::vector<int> busy_at_idle;
        void wait_events(std::vector<InputEvent>& out, int timeout_ms, int fd) override {
            if (timeout_ms < 0) busy_at_idle.push_back(busy_waits);
            FakeBackend::wait_events(out, timeout_ms, fd);
        }
    } be(kW, kH);
    be.script = {
        {{T::Press, p.x + 100, p.y + 20, 1}, {T::Release, p.x + 200, p.y + 150, 1}, {T::Leave}},  // zoom
        {{T::Press, p.x + 100, p.y + 20, 3}, {T::Leave}},                                         // unzoom
    };
    run_plot({sig}, o, be);
    ASSERT_EQ(be.busy_at_idle.size(), 3u);
    EXPECT_GT(be.busy_at_idle[1], be.busy_at_idle[0]) << "zooming in must compute";
    EXPECT_EQ(be.busy_at_idle[2], be.busy_at_idle[1]) << "unzooming must reuse, not compute";
    EXPECT_TRUE(same_plot_area(be.idle_frames[0], be.idle_frames[2], p));
}

TEST(Session, ModeChangeWhileZoomedRecomputesEveryLevel) {
    const auto v = sine(200'000);
    std::vector<std::complex<double>> c(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) c[i] = {v[i], 0.5 * v[(i * 7) % v.size()]};
    const Rect p = compute_layout(kW, kH).plot;
    FakeBackend be(kW, kH);
    InputEvent real{T::Key};
    real.key = '3';
    be.script = {
        {{T::Press, p.x + 100, p.y + 20, 1}, {T::Release, p.x + 200, p.y + 150, 1}, {T::Leave}},
        {real},
        {{T::Press, p.x + 100, p.y + 20, 3}, {T::Leave}},
    };
    run_plot({Signal(c)}, {}, be);
    ASSERT_EQ(be.idle_frames.size(), 4u);
    PlotOptions want;
    want.cmode = CMode::Real;
    EXPECT_TRUE(same_plot_area(be.idle_frames[3], headless({Signal(c)}, want, kW, kH), p))
        << "home level must show Real, not the cached Magnitude picture";
}

TEST(Session, ClickWhileRenderingStopsAndSpaceResumes) {
    // Large enough that an unoptimised build is still rendering when the click arrives.
    std::vector<float> v(60'000'000);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = static_cast<float>(i % 1000);
    const std::vector<Signal> sigs{Signal(v)};
    const Rect p = compute_layout(kW, kH).plot;
    FakeBackend be(kW, kH);
    be.while_busy = {{T::Press, p.x + 10, p.y + 10, 3}, {T::Leave}};  // right-click at home: just stops
    InputEvent space{T::Key};
    space.key = key::kSpace;
    be.script = {{space}};
    run_plot(sigs, {}, be);

    ASSERT_EQ(be.idle_frames.size(), 2u);
    const Framebuffer full = headless(sigs, {}, kW, kH);
    EXPECT_FALSE(same_plot_area(be.idle_frames[0], full, p)) << "stopped frame should be partial";
    EXPECT_TRUE(same_plot_area(be.idle_frames[1], full, p)) << "resumed frame should be complete";
}

TEST(Session, MenuAndHelpLeaveNoTraceWhenClosed) {
    const auto v = sine(50'000);
    const std::vector<Signal> sigs{Signal(v, 0.0, 1.0, "s")};
    const Rect p = compute_layout(kW, kH).plot;
    auto k = [](uint32_t sym) {
        InputEvent e{T::Key};
        e.key = sym;
        return e;
    };
    FakeBackend be(kW, kH);
    be.script = {
        {{T::Motion, p.x + 50, p.y + 50}, k('m')},  // menu at the pointer, over the plot
        {k(key::kDown), k(key::kRight)},           // open a submenu too
        {k(key::kEscape)},
        {k('?')},
        {k('a'), {T::Leave}},                       // any key closes help
    };
    run_plot(sigs, {}, be);
    ASSERT_EQ(be.idle_frames.size(), 6u);
    EXPECT_FALSE(same_plot_area(be.idle_frames[0], be.idle_frames[2], p)) << "menu is visible";
    EXPECT_FALSE(same_plot_area(be.idle_frames[0], be.idle_frames[4], p)) << "help is visible";
    EXPECT_TRUE(same_plot_area(be.idle_frames[0], be.idle_frames[5], p));
}

TEST(Session, IrKeyMatchesHeadlessIr) {
    std::vector<std::complex<float>> v(300'000);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = std::polar(1.0f + 0.1f * (i % 3), 0.001f * i);
    const std::vector<Signal> sigs{Signal(v)};
    InputEvent five{T::Key};
    five.key = '5';
    FakeBackend be(kW, kH);
    be.script = {{five}};
    run_plot(sigs, {}, be);
    ASSERT_EQ(be.idle_frames.size(), 2u);
    PlotOptions ir;
    ir.cmode = CMode::IR;
    EXPECT_TRUE(same_plot_area(be.idle_frames[1], headless(sigs, ir, kW, kH), compute_layout(kW, kH).plot));
}

std::vector<float> raster_data(std::size_t rows, std::size_t cols) {
    std::vector<float> v(rows * cols);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) v[r * cols + c] = std::sin(r * 0.05f) * std::cos(c * 0.07f) + r * 0.01f;
    return v;
}

TEST(SessionRaster, WindowMatchesHeadlessRaster) {
    const auto v = raster_data(300, 400);
    RasterOptions o;
    o.subsize = 400;
    FakeBackend be(kW, kH);
    run_session(std::make_unique<RasterContent>(Signal(v), o), be);
    ASSERT_EQ(be.idle_frames.size(), 1u);
    Screen s(std::make_unique<RasterContent>(Signal(v), o));
    Framebuffer want(kW, kH);
    render_headless(want, s);
    EXPECT_TRUE(same_plot_area(be.idle_frames[0], want, compute_layout(kW, kH, true).plot));
}

TEST(SessionRaster, CutIsIsolatedAndEscReturnsInstantly) {
    const auto v = raster_data(300, 400);
    RasterOptions o;
    o.subsize = 400;
    const Rect p = compute_layout(kW, kH, true).plot;
    auto key = [](uint32_t k) {
        InputEvent e{T::Key};
        e.key = k;
        return e;
    };
    class Counting : public FakeBackend {
    public:
        using FakeBackend::FakeBackend;
        std::vector<int> busy_at_idle;
        void wait_events(std::vector<InputEvent>& out, int timeout_ms, int fd) override {
            if (timeout_ms < 0) busy_at_idle.push_back(busy_waits);
            FakeBackend::wait_events(out, timeout_ms, fd);
        }
    } be(kW, kH);
    be.script = {
        {{T::Press, p.x + 40, p.y + 30, 1}, {T::Release, p.x + 200, p.y + 150, 1}, {T::Leave}},  // zoom raster
        {{T::Motion, p.x + 100, p.y + 60}, key('x'), {T::Leave}},                               // x-cut
        {{T::Press, p.x + 50, p.y + 20, 1}, {T::Release, p.x + 150, p.y + 120, 1}, key('3'), {T::Leave}},  // mess with the cut
        {key(0xff56), {T::Leave}},                                                               // PgDn
        {key(0xff55), {T::Leave}},                                                               // PgUp
        {key(key::kEscape)},                                                                     // leave the cut
    };
    run_session(std::make_unique<RasterContent>(Signal(v), o), be);
    ASSERT_EQ(be.idle_frames.size(), 7u);
    const auto& f = be.idle_frames;
    EXPECT_FALSE(same_plot_area(f[1], f[2], p)) << "the cut is showing";
    EXPECT_FALSE(same_plot_area(f[3], f[4], p)) << "PgDn shows the next row";
    EXPECT_TRUE(same_plot_area(f[3], f[5], p)) << "PgUp comes back, keeping the cut's zoom";
    EXPECT_TRUE(same_plot_area(f[1], f[6], p)) << "the raster is exactly as it was";
    EXPECT_EQ(be.busy_at_idle[6], be.busy_at_idle[5]) << "returning to the raster recomputes nothing";
}

TEST(Session, ResizeReflowsToTheNewSize) {
    const auto v = sine(10'000);
    FakeBackend be(kW, kH);
    InputEvent resize{T::Resize};
    resize.w = 700;
    resize.h = 400;
    be.script = {{resize}};
    run_plot({Signal(v)}, {}, be);
    ASSERT_EQ(be.idle_frames.size(), 2u);
    EXPECT_EQ(be.idle_frames[1].width(), 700);
    const Framebuffer want = headless({Signal(v)}, {}, 700, 400);
    const bool same = same_plot_area(be.idle_frames[1], want, compute_layout(700, 400).plot);
    if (!same && std::getenv("TEST_UNDECLARED_OUTPUTS_DIR")) {
        const std::string dir = std::getenv("TEST_UNDECLARED_OUTPUTS_DIR");
        write_png(be.idle_frames[1], dir + "/resize_actual.png");
        write_png(want, dir + "/resize_want.png");
    }
    EXPECT_TRUE(same);
}

TEST(Session, InterruptCheckEndsTheLoop) {
    const auto v = sine(1000);
    FakeBackend be(kW, kH);
    be.script = {{}, {}, {}};
    int calls = 0;
    run_plot({Signal(v)}, {}, be, [&] { return ++calls >= 2; });
    EXPECT_EQ(calls, 2);
}

TEST(JobRunner, NewJobCancelsTheRunningOne) {
    JobRunner jobs;
    std::atomic<bool> first_cancelled{false}, second_ran{false};
    std::atomic<bool> first_started{false};
    jobs.submit([&](const CancelToken& ct) {
        first_started = true;
        while (!ct.cancelled()) std::this_thread::yield();
        first_cancelled = true;
    });
    while (!first_started) std::this_thread::yield();
    jobs.submit([&](const CancelToken&) { second_ran = true; });
    jobs.wait_idle();
    EXPECT_TRUE(first_cancelled);
    EXPECT_TRUE(second_ran);

    pollfd fd{jobs.wake_fd(), POLLIN, 0};
    EXPECT_EQ(::poll(&fd, 1, 1000), 1);  // finishing a job wakes the UI
    jobs.drain();
    EXPECT_EQ(::poll(&fd, 1, 0), 0);
}

}  // namespace
}  // namespace ssp
