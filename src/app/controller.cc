#include "app/controller.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "app/menu.h"

namespace ssp {
namespace {

constexpr int kMinDragPx = 3;  // smaller boxes are clicks

// SigPlot's mode order: 1 Ma, 2 Ph, 3 Re, 4 Im, 5 IR, 6 Lo, 7 L2.
constexpr CMode kModeKeys[] = {CMode::Mag,  CMode::Phase, CMode::Real, CMode::Imag,
                               CMode::IR,   CMode::Log10, CMode::Log20};

/// Data coordinate of a pixel edge (not centre), for zoom boxes.
double edge_x(const Screen& s, int px) {
    const Range x = s.views.top().x;
    return x.lo + double(px - s.layout.plot.x) / s.layout.plot.w * (x.hi - x.lo);
}
double edge_y(const Screen& s, int py) {
    const Range y = s.yshown;
    return y.hi - double(py - s.layout.plot.y) / s.layout.plot.h * (y.hi - y.lo);
}

constexpr uint32_t kBackSpace = 0xff08, kEnter = 0xff0d, kKpEnter = 0xff8d;
constexpr std::size_t kMaxPromptChars = 60;

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
        s.ui.prompt.reset();
        View v = s.views.top();
        if (x) {
            v.x = r;
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
            s.content->cycle_style(t);
            // Dots need exact pixel rows, which only a reduce produces.
            return s.content->signal(t).style == Style::Dots ? kReduce : kRepaint;
        }
        return kNone;
    }
    if (!in_plot) return kNone;

    unsigned inv = rendering ? kCancel : kNone;  // any press in the plot stops rendering
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
    if (e.button != 1 || !s.ui.dragging) return kNone;
    s.ui.dragging = false;
    const Rect& p = s.layout.plot;
    const int xa = std::clamp(std::min(s.ui.x0, e.x), p.x, p.right());
    const int xb = std::clamp(std::max(s.ui.x0, e.x), p.x, p.right());
    const int ya = std::clamp(std::min(s.ui.y0, e.y), p.y, p.bottom());
    const int yb = std::clamp(std::max(s.ui.y0, e.y), p.y, p.bottom());
    if (xb - xa <= kMinDragPx || yb - ya <= kMinDragPx) return kOverlay;
    return action::zoom_to(s, View{{edge_x(s, xa), edge_x(s, xb)}, Range{edge_y(s, yb), edge_y(s, ya)}});
}

unsigned on_key(Screen& s, const InputEvent& e, bool rendering) {
    const uint32_t k = e.key;
    if (k == key::kEscape) {
        if (s.ui.dragging) {
            s.ui.dragging = false;
            return kOverlay;
        }
        return rendering ? kCancel : kNone;
    }
    if (k == key::kSpace) return kResume;
    if (k == key::kHome) return action::unzoom_all(s);
    if (k >= '1' && k <= '7') return action::set_mode(s, kModeKeys[k - '1']);
    if ((e.mods & mods::kCtrl) && (k == 's' || k == 'S')) return kSave;
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
    if (m == CMode::IR) {
        s.message = "IR (imag vs real) mode is not available yet";
        return kOverlay;
    }
    if (m == s.set.cmode) return kNone;
    s.set.cmode = m;
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
    const bool show = !s.content->signal(t).visible;
    s.content->set_visible(t, show);
    return show ? kReduce : kRepaint;  // hidden traces have no bins
}

unsigned toggle_index(Screen& s) {
    s.set.index = !s.set.index;  // different x coordinates: back to home
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

}  // namespace action
}  // namespace ssp
