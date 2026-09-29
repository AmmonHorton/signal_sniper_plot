#include "ssp/app/controller.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "ssp/app/menu.h"
#include "ssp/app/raster_content.h"
#include "ssp/render/colormap.h"

namespace ssp {
namespace {

constexpr int kMinDragPx = 3;  // smaller boxes are clicks

// SigPlot's mode order: 1 Ma, 2 Ph, 3 Re, 4 Im, 5 IR, 6 Lo, 7 L2.
constexpr CMode kModeKeys[] = {CMode::Mag,  CMode::Phase, CMode::Real, CMode::Imag,
                               CMode::IR,   CMode::Log10, CMode::Log20};

/// Data coordinate of a pixel edge (not centre), for zoom boxes.
double edge_x(const Screen& s, int px) {
    const Range x = s.xshown;
    return x.lo + double(px - s.layout.plot.x) / s.layout.plot.w * (x.hi - x.lo);
}
double edge_y(const Screen& s, int py) {
    const Range y = s.yshown;
    const double f = double(py - s.layout.plot.y) / s.layout.plot.h * (y.hi - y.lo);
    return s.ydown() ? y.lo + f : y.hi - f;
}

constexpr uint32_t kBackSpace = 0xff08, kEnter = 0xff0d, kKpEnter = 0xff8d;
constexpr uint32_t kPageUp = 0xff55, kPageDown = 0xff56;
constexpr double kZShift = 0.1;  // [ and ] slide the colour range by 10%
constexpr std::size_t kMaxPromptChars = 60;
constexpr double kWheelStep = 0.8;  // each wheel notch keeps 80% of the range (SigPlot: 20%)

/// Scale `r` by `f` about `anchor`, which stays at the same place on screen.
Range scale_about(Range r, double anchor, double f) {
    return {anchor - (anchor - r.lo) * f, anchor + (r.hi - anchor) * f};
}

/// Replace the adjustable top level with `v`, or push `v` as a new adjustable level.
unsigned adjust_view(Screen& s, View v) {
    v.adjustable = true;
    v.auto_x = false;
    v.cache.reset();
    if (s.views.top().adjustable) {
        s.views.top_mut() = v;
        return kReduce;
    }
    return action::zoom_to(s, v);
}

unsigned on_wheel(Screen& s, const InputEvent& e) {
    const bool in = e.button == 4;
    const bool do_x = !(e.mods & mods::kShift) || (e.mods & mods::kCtrl);
    const bool do_y = (e.mods & (mods::kShift | mods::kCtrl)) != 0;
    const double f = in ? kWheelStep : 1.0 / kWheelStep;
    View v = s.views.top();
    if (!in && !v.adjustable) return kNone;  // zooming out only undoes wheel zooms
    v.x = s.xshown;
    if (do_x) v.x = scale_about(s.xshown, s.px_to_x(e.x), f);
    if (do_y) v.y = scale_about(s.yshown, s.py_to_y(e.y), f);
    if (!in) {  // back past the level below: just return to it
        const View& below = s.views.parent();
        const Range bx = s.level_x(below);
        const bool x_out = !do_x || (v.x.lo <= bx.lo && v.x.hi >= bx.hi);
        const bool y_out = !do_y || !below.y || (v.y->lo <= below.y->lo && v.y->hi >= below.y->hi);
        if (x_out && y_out) return s.views.pop() ? kReduce : kNone;
    }
    return adjust_view(s, v);
}

/// Shift the view by (fx, fy) screen-widths/heights, keeping x on the data where possible.
unsigned pan_by(Screen& s, Range x0, Range y0, double fx, double fy) {
    View v = s.views.top();
    const double w = x0.hi - x0.lo, h = y0.hi - y0.lo;
    v.x = s.xshown;
    v.auto_x = false;
    if (s.set.cmode == CMode::IR) {  // no x extent to stay inside: pan freely
        if (fx != 0.0) v.x = {x0.lo + fx * w, x0.hi + fx * w};
    } else {
        const Range ext = s.home().x;
        if (fx != 0.0 && w < ext.hi - ext.lo) {  // x pans only within the data (SigPlot pan limits)
            v.x = {x0.lo + fx * w, x0.hi + fx * w};
            const double shift = std::max(ext.lo - v.x.lo, 0.0) - std::max(v.x.hi - ext.hi, 0.0);
            v.x = {v.x.lo + shift, v.x.hi + shift};
        }
    }
    if (fy != 0.0) v.y = Range{y0.lo + fy * h, y0.hi + fy * h};
    const View& cur = s.views.top();
    const bool same_y = v.y.has_value() == cur.y.has_value() && (!v.y || (v.y->lo == cur.y->lo && v.y->hi == cur.y->hi));
    const Range cx = s.level_x(cur);
    if (v.x.lo == cx.lo && v.x.hi == cx.hi && same_y) return kNone;
    return adjust_view(s, v);
}

/// "lo hi" (also "lo,hi" or "lo:hi") → range, or an error message.
bool parse_range(std::string text, Range& out) {
    std::replace_if(text.begin(), text.end(), [](char c) { return c == ',' || c == ':' || c == ';'; }, ' ');
    std::istringstream in(text);
    std::string rest;
    if (!(in >> out.lo >> out.hi) || (in >> rest)) return false;
    return std::isfinite(out.lo) && std::isfinite(out.hi) && out.hi > out.lo;
}

unsigned on_prompt_key(Screen& s, const InputEvent& e) {
    PromptState& p = *s.ui.prompt;
    if (e.key == key::kEscape) {
        s.ui.prompt.reset();
    } else if (e.key == kBackSpace) {
        if (!p.text.empty()) p.text.pop_back();
    } else if (e.key == kEnter || e.key == kKpEnter) {
        Range r;
        if (!parse_range(p.text, r)) {
            p.error = "enter two numbers, min then max";
            return kOverlay;
        }
        const bool x = p.kind == PromptState::Kind::XRange;
        if (p.kind == PromptState::Kind::ZRange) {
            s.ui.prompt.reset();
            s.set.zfixed = r;
            return kRepaint | kOverlay;
        }
        s.ui.prompt.reset();
        View v = s.views.top();
        if (x) {
            v.x = r;
            v.auto_x = false;
            v.y.reset();  // autoscale y over the new x range
        } else {
            v.y = r;
        }
        return action::zoom_to(s, v) | kOverlay;
    } else if (e.key >= 0x20 && e.key <= 0x7e && p.text.size() < kMaxPromptChars) {
        p.text += static_cast<char>(e.key);
        p.error.clear();
    }
    return kOverlay;
}

unsigned on_press(Screen& s, const InputEvent& e, bool rendering) {
    const bool in_plot = s.layout.plot.contains(e.x, e.y);
    if (e.button == 2) return ((rendering && in_plot) ? kCancel : kNone) | open_menu(s, e.x, e.y);
    for (std::size_t t = 0; t < s.legend_hits.size(); ++t) {
        if (!s.legend_hits[t].contains(e.x, e.y)) continue;
        if (e.button == 1) return action::toggle_trace(s, t);
        if (e.button == 3) {
            s.traces()->cycle_style(t);
            // Dots need exact pixel rows, which only a reduce produces.
            return s.traces()->signal(t).style == Style::Dots ? kReduce : kRepaint;
        }
        return kNone;
    }
    if (!in_plot) return kNone;

    unsigned inv = rendering ? kCancel : kNone;  // any press in the plot stops rendering
    if (e.button == 4 || e.button == 5) return inv | on_wheel(s, e);
    if (e.button == 1 && (e.mods & mods::kShift)) {
        s.ui.panning = true;
        s.ui.x0 = e.x;
        s.ui.y0 = e.y;
        s.ui.pan_x = s.xshown;
        s.ui.pan_y = s.yshown;
        return inv | kOverlay;
    }
    if (e.button == 1) {
        s.ui.dragging = true;
        s.ui.x0 = e.x;
        s.ui.y0 = e.y;
        inv |= kOverlay;
    } else if (e.button == 3) {
        if (s.views.pop()) inv |= kReduce;
    }
    return inv;
}

unsigned on_release(Screen& s, const InputEvent& e) {
    if (e.button == 1 && s.ui.panning) {
        s.ui.panning = false;
        return kOverlay;
    }
    if (e.button != 1 || !s.ui.dragging) return kNone;
    s.ui.dragging = false;
    const Rect& p = s.layout.plot;
    const int xa = std::clamp(std::min(s.ui.x0, e.x), p.x, p.right());
    const int xb = std::clamp(std::max(s.ui.x0, e.x), p.x, p.right());
    const int ya = std::clamp(std::min(s.ui.y0, e.y), p.y, p.bottom());
    const int yb = std::clamp(std::max(s.ui.y0, e.y), p.y, p.bottom());
    if (std::abs(e.x - s.ui.x0) <= kMinDragPx && std::abs(e.y - s.ui.y0) <= kMinDragPx) {
        s.ui.marker = {{s.px_to_x(e.x), s.py_to_y(e.y)}};  // a click: set the marker
        return kOverlay;
    }
    if (xb - xa <= kMinDragPx || yb - ya <= kMinDragPx) return kOverlay;
    const double y0 = edge_y(s, ya), y1 = edge_y(s, yb);
    return action::zoom_to(s, View{{edge_x(s, xa), edge_x(s, xb)}, Range{std::min(y0, y1), std::max(y0, y1)}});
}

unsigned on_key(Screen& s, const InputEvent& e, bool rendering) {
    const uint32_t k = e.key;
    if (k == key::kEscape) {
        if (s.ui.dragging) {
            s.ui.dragging = false;
            return kOverlay;
        }
        s.ui.panning = false;
        if (rendering) return kCancel;
        return s.cut ? kPop : kNone;  // leave the cut, back to the raster
    }
    if (k == key::kSpace) return kResume;
    if (k == key::kHome) return action::unzoom_all(s);
    if (k == key::kLeft || k == key::kRight) return pan_by(s, s.xshown, s.yshown, k == key::kLeft ? -0.5 : 0.5, 0.0);
    if (k == key::kUp || k == key::kDown) {  // "up" shows what is above on screen
        const double up = s.ydown() ? -0.5 : 0.5;
        return pan_by(s, s.xshown, s.yshown, 0.0, k == key::kUp ? up : -up);
    }
    if (k >= '1' && k <= '7') return action::set_mode(s, kModeKeys[k - '1']);
    if ((e.mods & mods::kCtrl) && (k == 's' || k == 'S')) return kSave;
    if ((k == kPageUp || k == kPageDown) && s.cut) {
        s.cut_step = k == kPageDown ? 1 : -1;
        return kStep;
    }
    if (s.content->is_raster()) {
        switch (k) {
            case 'x':
            case 'y':
                if (!s.ui.inside) return kNone;
                return action::request_cut(s, k == 'x', s.ui.mx, s.ui.my);
            case 'c':
                return action::cycle_colormap(s, 1);
            case 'C':
                return action::cycle_colormap(s, -1);
            case '[':
                return action::shift_z(s, -kZShift);
            case ']':
                return action::shift_z(s, kZShift);
            default:
                break;
        }
    }
    switch (k) {
        case 'i':
            return action::toggle_index(s);
        case 'm':
            return s.ui.inside ? open_menu(s, s.ui.mx, s.ui.my)
                               : open_menu(s, s.layout.plot.x + 20, s.layout.plot.y + 20);
        case '?':
            s.ui.help = true;
            return kOverlay;
        case 'g':
            s.set.grid = !s.set.grid;
            return kRepaint;
        case 'k':
            s.ui.marker.reset();
            return kOverlay;
        case 'a':
            s.set.absc = s.set.absc == Absc::X ? Absc::Index
                       : s.set.absc == Absc::Index ? Absc::Inverse : Absc::X;
            return kOverlay;
        case 'l':
            s.set.legend = !s.set.legend;
            return kRepaint;
        case 'q':
            return kQuit;
        default:
            return kNone;
    }
}

}  // namespace

unsigned handle_event(Screen& s, const InputEvent& e, bool rendering) {
    using T = InputEvent::Type;
    if (e.type == T::Close) return kQuit;
    if (e.type == T::Motion || e.type == T::Press || e.type == T::Release) {
        s.ui.mx = e.x;
        s.ui.my = e.y;
        s.ui.inside = true;
    }
    // Popups take input first: help closes on any key or click, the menu takes everything.
    if (s.ui.help && (e.type == T::Key || e.type == T::Press)) {
        s.ui.help = false;
        return kOverlay;
    }
    if (s.ui.menu && e.type != T::Leave) return handle_menu_event(s, e) | kOverlay;
    if (s.ui.prompt && e.type == T::Key) return on_prompt_key(s, e);

    switch (e.type) {
        case T::Motion:
            if (s.ui.panning) {  // grab-and-drag: the point under the pointer follows it
                const double fx = -double(e.x - s.ui.x0) / s.layout.plot.w;
                const double fy = (s.ydown() ? -1.0 : 1.0) * double(e.y - s.ui.y0) / s.layout.plot.h;
                return pan_by(s, s.ui.pan_x, s.ui.pan_y, fx, fy) | kOverlay;
            }
            return kOverlay;
        case T::Leave:
            s.ui.inside = false;
            return kOverlay;
        case T::Press:
            s.message.clear();
            return on_press(s, e, rendering) | kOverlay;
        case T::Release:
            return on_release(s, e);
        case T::Key:
            s.message.clear();
            return on_key(s, e, rendering) | kOverlay;
        default:
            return kNone;
    }
}

namespace action {

unsigned set_mode(Screen& s, CMode m) {
    if (m == s.set.cmode) return kNone;
    if (!s.content->supports(m)) {
        s.message = std::string(cmode_name(m)) + " is not available here";
        return kOverlay;
    }
    const bool was_ir = s.set.cmode == CMode::IR;
    // Entering IR plots the samples of the time window in view now.
    if (m == CMode::IR) s.ir = s.traces()->samples_in(s.xshown, s.set.index);
    s.set.cmode = m;
    s.ui.marker.reset();  // its y was in the old units
    if (m == CMode::IR || was_ir) {  // different axes altogether: start over at home
        s.go_home();
        return kReduce;
    }
    // Same x at every level, y re-autoscaled (SigPlot changemode); a fixed yrange only
    // applies to the mode it was given for.
    s.views.clear_y();
    s.views.clear_cache();  // every level was computed for the old mode
    s.views.home().y = s.home().y;
    return kReduce;
}

unsigned set_phunits(Screen& s, PhaseUnits u) {
    if (u == s.set.phunits) return kNone;
    s.set.phunits = u;
    if (s.set.cmode != CMode::Phase) return kNone;
    s.views.clear_y();  // phase values changed scale
    s.views.clear_cache();
    return kReduce;
}

unsigned toggle_trace(Screen& s, std::size_t t) {
    const bool show = !s.traces()->signal(t).visible;
    s.traces()->set_visible(t, show);
    return show ? kReduce : kRepaint;  // hidden traces have no bins
}

unsigned set_style(Screen& s, Style st) {
    TraceContent* tc = s.traces();
    if (!tc) return kNone;
    tc->set_style(st);
    return st == Style::Dots ? kReduce : kRepaint;  // dots need exact pixel rows
}

unsigned toggle_index(Screen& s) {
    s.set.index = !s.set.index;  // different x coordinates: back to home
    s.ui.marker.reset();
    s.go_home();
    return kReduce;
}

unsigned unzoom_all(Screen& s) {
    if (s.views.level() == 0) return kNone;
    s.views.pop_all();
    return kReduce;
}

unsigned autoscale_y(Screen& s) {
    if (!s.views.top().y) return kNone;
    View v = s.views.top();
    v.y.reset();
    return zoom_to(s, v);
}

unsigned zoom_to(Screen& s, const View& v) {
    View level = v;
    level.cache.reset();  // a new level starts with nothing computed
    if (!s.views.push(level)) {
        s.message = "zoom limit reached (10 levels)";
        return kOverlay;
    }
    return kReduce;
}

unsigned open_prompt(Screen& s, PromptState::Kind kind) {
    s.ui.prompt = PromptState{kind, {}, {}};
    return kOverlay;
}

unsigned request_cut(Screen& s, bool xcut, int px, int py) {
    const RasterContent* rc = s.raster();
    if (!rc) return kNone;
    const auto cell = xcut ? rc->row_at(s.py_to_y(py), s.set.index) : rc->col_at(s.px_to_x(px), s.set.index);
    if (!s.layout.plot.contains(px, py) || !cell) {
        s.message = "point at the raster first";
        return kOverlay;
    }
    s.cut_request = {{xcut, *cell}};
    return kCut;
}

unsigned cycle_colormap(Screen& s, int step) {
    const int n = kNumColormaps;
    s.set.cmap = static_cast<Colormap>((static_cast<int>(s.set.cmap) + step + n) % n);
    s.message = std::string("colormap: ") + colormap_name(s.set.cmap);
    return kRepaint;
}

unsigned set_reduce(Screen& s, Reduce r) {
    if (r == s.set.reduce) return kNone;
    s.set.reduce = r;
    return kReduce;
}

unsigned shift_z(Screen& s, double f) {
    const Range z = s.zshown;
    const double d = f * (z.hi - z.lo);
    s.set.zfixed = Range{z.lo + d, z.hi + d};
    return kRepaint;
}

}  // namespace action
}  // namespace ssp
