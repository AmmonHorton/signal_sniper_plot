#include "app/controller.h"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "app/compose.h"
#include "app/menu.h"

namespace ssp {
namespace {

using T = InputEvent::Type;

InputEvent press(int x, int y, int button = 1) { return {T::Press, x, y, button}; }
InputEvent release(int x, int y, int button = 1) { return {T::Release, x, y, button}; }
InputEvent key(uint32_t k, unsigned mods = 0) {
    InputEvent e{T::Key};
    e.key = k;
    e.mods = mods;
    return e;
}

/// A screen that has been rendered once, so layout, yshown and legend hits are real.
class ControllerTest : public ::testing::Test {
protected:
    void SetUp() override { make({Signal(a_, 0.0, 1.0, "a"), Signal(b_, 0.0, 1.0, "b")}); }

    void make(std::vector<Signal> sigs, PlotOptions o = {}) {
        screen_ = std::make_unique<Screen>(std::make_unique<TraceContent>(std::move(sigs), o));
        render_headless(fb_, *screen_);
    }
    Screen& s() { return *screen_; }
    const Rect& plot() { return screen_->layout.plot; }

    std::vector<std::complex<float>> a_ = std::vector<std::complex<float>>(1000, {1.0f, 2.0f});
    std::vector<std::complex<float>> b_ = std::vector<std::complex<float>>(1000, {-1.0f, 0.5f});
    Framebuffer fb_{500, 300};
    std::unique_ptr<Screen> screen_;
};

TEST_F(ControllerTest, DragZoomsToTheBoxAndRightClickUnzooms) {
    const Range x0 = s().views.top().x;
    const int left = plot().x + plot().w / 4, right = plot().x + plot().w / 2;
    EXPECT_EQ(handle_event(s(), press(left, plot().y + 10), false) & kReduce, 0u);
    EXPECT_TRUE(handle_event(s(), release(right, plot().y + 60), false) & kReduce);
    ASSERT_EQ(s().views.level(), 1u);
    const View v = s().views.top();
    const double w = x0.hi - x0.lo;
    EXPECT_NEAR(v.x.lo, x0.lo + 0.25 * w, w / plot().w);
    EXPECT_NEAR(v.x.hi, x0.lo + 0.5 * w, w / plot().w);
    ASSERT_TRUE(v.y.has_value());
    EXPECT_GT(v.y->hi, v.y->lo);

    EXPECT_TRUE(handle_event(s(), press(left, plot().y + 10, 3), false) & kReduce);
    EXPECT_EQ(s().views.level(), 0u);
    EXPECT_EQ(handle_event(s(), press(left, plot().y + 10, 3), false) & kReduce, 0u);  // already home
}

TEST_F(ControllerTest, TinyDragIsNotAZoom) {
    handle_event(s(), press(plot().x + 50, plot().y + 50), false);
    EXPECT_EQ(handle_event(s(), release(plot().x + 52, plot().y + 90), false) & kReduce, 0u);
    EXPECT_EQ(s().views.level(), 0u);
}

TEST_F(ControllerTest, AnyPressInThePlotStopsRendering) {
    EXPECT_TRUE(handle_event(s(), press(plot().x + 5, plot().y + 5, 2), true) & kCancel);
    EXPECT_FALSE(handle_event(s(), press(plot().x + 5, plot().y + 5, 2), false) & kCancel);
    EXPECT_TRUE(handle_event(s(), key(key::kEscape), true) & kCancel);
}

TEST_F(ControllerTest, EscapeAbandonsADrag) {
    handle_event(s(), press(plot().x + 50, plot().y + 50), false);
    ASSERT_TRUE(s().ui.dragging);
    handle_event(s(), key(key::kEscape), false);
    EXPECT_FALSE(s().ui.dragging);
    EXPECT_EQ(handle_event(s(), release(plot().x + 200, plot().y + 150), false) & kReduce, 0u);
}

TEST_F(ControllerTest, ModeKeysKeepXAndRescaleYAtEveryLevel) {
    handle_event(s(), press(plot().x + 20, plot().y + 20), false);
    handle_event(s(), release(plot().x + 200, plot().y + 150), false);
    const Range zoomed_x = s().views.top().x;
    ASSERT_TRUE(s().views.top().y);

    EXPECT_TRUE(handle_event(s(), key('3'), false) & kReduce);
    EXPECT_EQ(s().set.cmode, CMode::Real);
    EXPECT_EQ(s().views.level(), 1u);
    EXPECT_EQ(s().views.top().x.lo, zoomed_x.lo);
    EXPECT_FALSE(s().views.top().y.has_value());

    EXPECT_EQ(handle_event(s(), key('3'), false) & kReduce, 0u);  // no change
}

TEST_F(ControllerTest, IrUsesTheTimeWindowInViewAndResetsLevels) {
    action::zoom_to(s(), View{{100.0, 199.5}, std::nullopt});
    s().xshown = {100.0, 199.5};
    EXPECT_TRUE(handle_event(s(), key('5'), false) & kReduce);
    EXPECT_EQ(s().set.cmode, CMode::IR);
    EXPECT_EQ(s().views.level(), 0u);
    EXPECT_TRUE(s().views.top().auto_x);
    ASSERT_EQ(s().ir.size(), 2u);
    EXPECT_EQ(s().ir[0], (SampleRange{100, 200}));  // samples whose x is in [100, 199.5]

    // Box zoom inside IR gives explicit axes; leaving IR returns to the time home.
    s().xshown = {-2.0, 2.0};
    s().yshown = {-2.0, 2.0};
    handle_event(s(), press(plot().x + 20, plot().y + 20), false);
    handle_event(s(), release(plot().x + 200, plot().y + 150), false);
    EXPECT_FALSE(s().views.top().auto_x);
    EXPECT_EQ(s().views.level(), 1u);
    EXPECT_TRUE(handle_event(s(), key('1'), false) & kReduce);
    EXPECT_EQ(s().views.level(), 0u);
    EXPECT_FALSE(s().views.top().auto_x);
    EXPECT_EQ(s().views.top().x.hi, 999.0);
}

TEST_F(ControllerTest, LegendTogglesVisibilityAndStyle) {
    ASSERT_EQ(s().legend_hits.size(), 2u);
    const Rect r = s().legend_hits[1];
    EXPECT_TRUE(handle_event(s(), press(r.x + 2, r.y + 2), false) & kRepaint);  // hide: repaint
    EXPECT_FALSE(s().traces()->signal(1).visible);
    EXPECT_TRUE(handle_event(s(), press(r.x + 2, r.y + 2), false) & kReduce);   // show: needs bins
    EXPECT_TRUE(s().traces()->signal(1).visible);
    EXPECT_TRUE(handle_event(s(), press(r.x + 2, r.y + 2, 3), false) & kReduce);  // → dots
    EXPECT_EQ(s().traces()->signal(1).style, Style::Dots);
    EXPECT_TRUE(handle_event(s(), press(r.x + 2, r.y + 2, 3), false) & kRepaint);  // → both
    EXPECT_EQ(s().traces()->signal(1).style, Style::LinesDots);
    EXPECT_EQ(s().views.level(), 0u);  // legend clicks never start a zoom
}

TEST_F(ControllerTest, IndexKeySwitchesAxesAndGoesHome) {
    make({Signal(a_, 100.0, 0.5, "a")});
    handle_event(s(), press(plot().x + 20, plot().y + 20), false);
    handle_event(s(), release(plot().x + 200, plot().y + 150), false);
    EXPECT_TRUE(handle_event(s(), key('i'), false) & kReduce);
    EXPECT_TRUE(s().set.index);
    EXPECT_EQ(s().views.level(), 0u);
    EXPECT_EQ(s().views.top().x.lo, 0.0);
    EXPECT_EQ(s().views.top().x.hi, 999.0);
}

TEST_F(ControllerTest, HomeUnzoomsAllLevels) {
    for (int i = 0; i < 3; ++i) {
        handle_event(s(), press(plot().x + 20, plot().y + 20), false);
        handle_event(s(), release(plot().x + 200, plot().y + 150), false);
    }
    ASSERT_EQ(s().views.level(), 3u);
    EXPECT_TRUE(handle_event(s(), key(key::kHome), false) & kReduce);
    EXPECT_EQ(s().views.level(), 0u);
}

TEST_F(ControllerTest, FixedYRangeOnlyAppliesToItsMode) {
    PlotOptions o;
    o.yrange = Range{-5, 5};
    make({Signal(a_)}, o);
    EXPECT_EQ(s().views.top().y->lo, -5.0);
    handle_event(s(), key('3'), false);
    EXPECT_FALSE(s().views.top().y.has_value());
    handle_event(s(), key('1'), false);
    ASSERT_TRUE(s().views.top().y.has_value());
    EXPECT_EQ(s().views.top().y->lo, -5.0);
}

TEST_F(ControllerTest, OtherKeys) {
    EXPECT_TRUE(handle_event(s(), key('s', mods::kCtrl), false) & kSave);
    EXPECT_TRUE(handle_event(s(), key('q'), false) & kQuit);
    EXPECT_TRUE(handle_event(s(), key(key::kSpace), false) & kResume);
    EXPECT_TRUE(handle_event(s(), key('g'), false) & kRepaint);
    EXPECT_FALSE(s().set.grid);
    EXPECT_TRUE(handle_event(s(), InputEvent{T::Close}, false) & kQuit);
}

InputEvent motion(int x, int y) { return {T::Motion, x, y}; }

TEST_F(ControllerTest, MenuByKeyboardSetsMode) {
    handle_event(s(), key('m'), false);
    ASSERT_TRUE(s().ui.menu.has_value());
    handle_event(s(), key(key::kDown), false);   // Mode
    handle_event(s(), key(key::kRight), false);  // → Magnitude
    handle_event(s(), key(key::kDown), false);   // Phase
    handle_event(s(), key(key::kDown), false);   // Real
    EXPECT_TRUE(handle_event(s(), key(0xff0d), false) & kReduce);
    EXPECT_EQ(s().set.cmode, CMode::Real);
    EXPECT_FALSE(s().ui.menu.has_value());
}

TEST_F(ControllerTest, MenuEatsInputUntilClosed) {
    const int x = plot().x + 30, y = plot().y + 30;
    handle_event(s(), press(x, y, 2), false);
    ASSERT_TRUE(s().ui.menu.has_value());
    handle_event(s(), motion(x + 20, y + 13), false);  // over the first item ("Mode")
    EXPECT_EQ(s().ui.menu->sel.size(), 2u) << "hovering a submenu item opens it";
    EXPECT_EQ(handle_event(s(), key('q'), false) & kQuit, 0u) << "keys go to the menu";
    // A click outside closes the menu and is not the start of a zoom.
    handle_event(s(), press(plot().right() - 5, plot().bottom() - 5), false);
    EXPECT_FALSE(s().ui.menu.has_value());
    EXPECT_FALSE(s().ui.dragging);
}

TEST_F(ControllerTest, MiddleClickWhileRenderingStopsAndOpensMenu) {
    const unsigned inv = handle_event(s(), press(plot().x + 30, plot().y + 30, 2), true);
    EXPECT_TRUE(inv & kCancel);
    EXPECT_TRUE(s().ui.menu.has_value());
}

void type(Screen& s, const std::string& text) {
    for (char c : text) handle_event(s, key(static_cast<uint8_t>(c)), false);
}

TEST_F(ControllerTest, RangePromptPushesAZoomLevel) {
    action::open_prompt(s(), PromptState::Kind::XRange);
    type(s(), "10, 20.5");
    EXPECT_TRUE(handle_event(s(), key(0xff0d), false) & kReduce);
    EXPECT_FALSE(s().ui.prompt.has_value());
    ASSERT_EQ(s().views.level(), 1u);
    EXPECT_EQ(s().views.top().x.lo, 10.0);
    EXPECT_EQ(s().views.top().x.hi, 20.5);
    EXPECT_FALSE(s().views.top().y.has_value()) << "y autoscales over the new x range";

    action::open_prompt(s(), PromptState::Kind::YRange);
    type(s(), "-1 1");
    handle_event(s(), key(0xff0d), false);
    ASSERT_EQ(s().views.level(), 2u);
    EXPECT_EQ(s().views.top().x.lo, 10.0);
    EXPECT_EQ(s().views.top().y->hi, 1.0);
}

TEST_F(ControllerTest, RangePromptRejectsBadInputAndEscCancels) {
    action::open_prompt(s(), PromptState::Kind::XRange);
    type(s(), "5 1");
    EXPECT_EQ(handle_event(s(), key(0xff0d), false) & kReduce, 0u);
    ASSERT_TRUE(s().ui.prompt.has_value());
    EXPECT_FALSE(s().ui.prompt->error.empty());
    handle_event(s(), key(0xff08), false);  // backspace
    EXPECT_EQ(s().ui.prompt->text, "5 ");
    handle_event(s(), key(key::kEscape), false);
    EXPECT_FALSE(s().ui.prompt.has_value());
    EXPECT_EQ(s().views.level(), 0u);
}

TEST_F(ControllerTest, HelpClosesOnTheNextKeyWhichIsConsumed) {
    handle_event(s(), key('?'), false);
    EXPECT_TRUE(s().ui.help);
    EXPECT_EQ(handle_event(s(), key('q'), false) & kQuit, 0u);
    EXPECT_FALSE(s().ui.help);
}

TEST_F(ControllerTest, PhaseUnitsOnlyRecomputeInPhaseMode) {
    EXPECT_EQ(action::set_phunits(s(), PhaseUnits::Degrees) & kReduce, 0u);
    action::set_mode(s(), CMode::Phase);
    EXPECT_TRUE(action::set_phunits(s(), PhaseUnits::Cycles) & kReduce);
}

InputEvent wheel(int x, int y, int button, unsigned mods = 0) {
    InputEvent e{T::Press, x, y, button};
    e.mods = mods;
    return e;
}

TEST_F(ControllerTest, ClickSetsMarkerAndKClearsIt) {
    const int x = plot().x + 100, y = plot().y + 50;
    handle_event(s(), press(x, y), false);
    EXPECT_EQ(handle_event(s(), release(x + 1, y), false) & kReduce, 0u);
    ASSERT_TRUE(s().ui.marker.has_value());
    EXPECT_NEAR(s().ui.marker->first, s().px_to_x(x + 1), 1e-12);
    EXPECT_EQ(s().views.level(), 0u);
    handle_event(s(), key('k'), false);
    EXPECT_FALSE(s().ui.marker.has_value());
    handle_event(s(), press(x, y), false);
    handle_event(s(), release(x, y), false);
    action::set_mode(s(), CMode::Real);
    EXPECT_FALSE(s().ui.marker.has_value()) << "mode change invalidates the marker's y";
}

TEST_F(ControllerTest, WheelZoomsAboutThePointerInOneLevel) {
    const int x = plot().x + plot().w / 4, y = plot().y + 40;
    const Range home = s().views.top().x;
    const double anchor = s().px_to_x(x);
    EXPECT_TRUE(handle_event(s(), wheel(x, y, 4), false) & kReduce);
    ASSERT_EQ(s().views.level(), 1u);
    const Range z = s().views.top().x;
    EXPECT_NEAR(z.hi - z.lo, 0.8 * (home.hi - home.lo), 1e-9);
    EXPECT_NEAR((anchor - z.lo) / (z.hi - z.lo), (anchor - home.lo) / (home.hi - home.lo), 1e-9)
        << "the point under the pointer stays put";

    s().xshown = s().views.top().x;  // what compose would show
    handle_event(s(), wheel(x, y, 4), false);
    EXPECT_EQ(s().views.level(), 1u) << "further wheel steps adjust the same level";
    for (int i = 0; i < 4; ++i) {
        s().xshown = s().views.top().x;
        handle_event(s(), wheel(x, y, 5), false);
    }
    EXPECT_EQ(s().views.level(), 0u) << "wheeling out returns to the level below";
    EXPECT_EQ(handle_event(s(), wheel(x, y, 5), false) & kReduce, 0u) << "never out past home";
}

TEST_F(ControllerTest, ShiftWheelZoomsYOnly) {
    const Range x0 = s().views.top().x;
    handle_event(s(), wheel(plot().x + 50, plot().y + 50, 4, mods::kShift), false);
    EXPECT_EQ(s().views.top().x.lo, x0.lo);
    ASSERT_TRUE(s().views.top().y.has_value());
    EXPECT_NEAR(s().views.top().y->hi - s().views.top().y->lo, 0.8 * (s().yshown.hi - s().yshown.lo), 1e-9);
}

TEST_F(ControllerTest, ArrowsPanWithinTheData) {
    EXPECT_EQ(handle_event(s(), key(key::kRight), false) & kReduce, 0u) << "home already shows all x";
    EXPECT_EQ(s().views.level(), 0u);
    action::zoom_to(s(), View{{100.0, 200.0}, std::nullopt});
    s().xshown = {100.0, 200.0};
    EXPECT_TRUE(handle_event(s(), key(key::kRight), false) & kReduce);
    EXPECT_EQ(s().views.top().x.lo, 150.0);
    EXPECT_EQ(s().views.level(), 2u);
    s().xshown = {900.0, 1000.0};
    handle_event(s(), key(key::kRight), false);
    EXPECT_EQ(s().views.top().x.hi, 999.0) << "clamped to the end of the data";
    EXPECT_EQ(s().views.level(), 2u) << "panning adjusts one level";
}

TEST_F(ControllerTest, ShiftDragPans) {
    action::zoom_to(s(), View{{100.0, 200.0}, Range{-1.0, 1.0}});
    s().xshown = {100.0, 200.0};
    s().yshown = {-1.0, 1.0};
    const int x = plot().x + 100, y = plot().y + 50;
    InputEvent p = press(x, y);
    p.mods = mods::kShift;
    handle_event(s(), p, false);
    EXPECT_FALSE(s().ui.dragging);
    EXPECT_TRUE(handle_event(s(), motion(x + plot().w / 10, y), false) & kReduce);
    const double moved = double(plot().w / 10) / plot().w * 100.0;
    EXPECT_NEAR(s().views.top().x.lo, 100.0 - moved, 1e-9) << "data follows the pointer";
    handle_event(s(), release(x + plot().w / 10, y), false);
    EXPECT_FALSE(s().ui.panning);
    EXPECT_EQ(s().views.level(), 2u);
}

TEST_F(ControllerTest, AKeyCyclesReadout) {
    EXPECT_EQ(s().set.absc, Absc::X);
    handle_event(s(), key('a'), false);
    EXPECT_EQ(s().set.absc, Absc::Index);
    handle_event(s(), key('a'), false);
    EXPECT_EQ(s().set.absc, Absc::Inverse);
    handle_event(s(), key('a'), false);
    EXPECT_EQ(s().set.absc, Absc::X);
}

TEST_F(ControllerTest, LegendKeyShowsASingleUnnamedTrace) {
    make({Signal(a_)});  // one unnamed trace: no legend at first
    EXPECT_FALSE(s().set.legend);
    EXPECT_TRUE(s().legend_hits.empty());
    EXPECT_TRUE(handle_event(s(), key('l'), false) & kRepaint);
    render_headless(fb_, s());
    ASSERT_EQ(s().legend_hits.size(), 1u) << "'l' shows it anyway";
    const Rect r = s().legend_hits[0];
    EXPECT_TRUE(handle_event(s(), press(r.x + 2, r.y + 2, 3), false) & kReduce);
    EXPECT_EQ(s().traces()->signal(0).style, Style::Dots);
}

TEST_F(ControllerTest, StyleMenuSetsEveryTrace) {
    handle_event(s(), key('m'), false);
    for (int i = 0; i < 4; ++i) handle_event(s(), key(key::kDown), false);  // Mode, Phase units, Scaling, Style
    handle_event(s(), key(key::kRight), false);                            // → Lines
    handle_event(s(), key(key::kDown), false);                             // Dots
    EXPECT_TRUE(handle_event(s(), key(0xff0d), false) & kReduce);
    EXPECT_EQ(s().traces()->signal(0).style, Style::Dots);
    EXPECT_EQ(s().traces()->signal(1).style, Style::Dots);
    EXPECT_TRUE(action::set_style(s(), Style::Lines) & kRepaint);
    EXPECT_EQ(s().traces()->signal(1).style, Style::Lines);
}

}  // namespace
}  // namespace ssp
