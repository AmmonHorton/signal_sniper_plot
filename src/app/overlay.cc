#include "app/overlay.h"

#include <algorithm>
#include <cmath>
#include <string>

#include "app/menu.h"

namespace ssp {
namespace {

constexpr int kCW = Framebuffer::kCharW;
constexpr int kCH = Framebuffer::kCharH;
constexpr int kCrossDot = 2;

// SigPlot readout abbreviations.
const char* cmode_short(CMode m) {
    switch (m) {
        case CMode::Mag:   return "Ma";
        case CMode::Phase: return "Ph";
        case CMode::Real:  return "Re";
        case CMode::Imag:  return "Im";
        case CMode::IR:    return "IR";
        case CMode::Log10: return "Lo";
        case CMode::Log20: return "L2";
        case CMode::Auto:  break;
    }
    return "??";
}

// Keypress help, SigPlot-style table. Keep in sync with controller.cc and menu.cc.
const char* const kHelp[] = {
    "Mouse",
    "  left drag        zoom to box",
    "  right click      unzoom one level",
    "  middle click     menu",
    "  legend click     show/hide trace (right click: lines/dots/both)",
    "  left click       set marker (readout shows dx, dy from it)",
    "  shift+left drag  pan",
    "  wheel            zoom x (shift: y, ctrl: both)",
    "  any click        stop a render in progress",
    "",
    "Keys",
    "  1-7              Mag Phase Real Imag ImagVsReal 10log 20log",
    "                   (5 = Imag vs Real of the time window in view)",
    "  Home             unzoom all",
    "  arrows           pan half a screen",
    "  a                readout: x / sample index / 1/x",
    "  k                clear marker",
    "  Space / Esc      resume / stop rendering",
    "  i                index x-axis",
    "  g  l             grid, legend",
    "                   (menu > Style: lines / dots for every trace)",
    "  m                menu (also: x/y range, phase units, traces)",
    "  Ctrl-S           save PNG",
    "  q                quit",
    "  ?                this help (any key closes)",
    "",
    "Raster",
    "  x / y            cut the row / column under the pointer (within the zoom box)",
    "  Esc              leave a cut, back to the raster",
    "  PgUp / PgDn      previous / next row or column, in a cut",
    "  c / C            next / previous colormap",
    "  [ / ]            slide the colour range down / up",
};

Rect draw_help(Framebuffer& fb, const Theme& th) {
    std::size_t chars = 0;
    for (const char* line : kHelp) chars = std::max(chars, std::string(line).size());
    constexpr int kLineH = kCH + 2, kPad = 8;
    const int n = static_cast<int>(sizeof kHelp / sizeof kHelp[0]);
    const int w = static_cast<int>(chars) * kCW + 2 * kPad;
    const int h = n * kLineH + 2 * kPad;
    const Rect box{std::max(0, (fb.width() - w) / 2), std::max(0, (fb.height() - h) / 2), w, h};
    fb.fill_rect(box, th.bg);
    fb.rect_outline(box, th.fg);
    for (int i = 0; i < n; ++i) fb.text(box.x + kPad, box.y + kPad + i * kLineH, kHelp[i], th.fg);
    return box.intersect(fb.bounds());
}

/// Readout text for an x (or dx) value in the chosen abscissa mode.
std::string x_text(const Screen& s, double x, bool delta) {
    const Range axis = s.content->xaxis();  // {xstart, xdelta} of the first trace / columns
    if (s.set.cmode == CMode::IR) return format_g(x, 9);  // x is the real part, not time
    switch (s.set.absc) {
        case Absc::Index: {
            // Sample index of the first trace (SigPlot uses the first layer too).
            const double i = s.set.index ? x : (x - (delta ? 0.0 : axis.lo)) / axis.hi;
            return std::to_string(std::llround(i));
        }
        case Absc::Inverse:
            return format_g(x == 0.0 ? 0.0 : 1.0 / x, 9);
        case Absc::X:
            break;
    }
    return s.set.index ? std::to_string(std::llround(x)) : format_g(x, 9);
}

std::string y_text(const Screen& s, double y) {
    return format_g(s.set.absc == Absc::Inverse && y != 0.0 ? 1.0 / y : y, 9);
}

const char* absc_label(Absc a) {
    switch (a) {
        case Absc::Index:   return "(indx)";
        case Absc::Inverse: return "(1/ab)";
        case Absc::X:       break;
    }
    return "(absc)";
}

std::string percent(double f) {
    return std::to_string(static_cast<int>(std::floor(100.0 * std::clamp(f, 0.0, 1.0)))) + "%";
}

}  // namespace

std::vector<Rect> draw_overlay(Framebuffer& fb, const Screen& s, const Theme& th, bool rendering,
                               double progress) {
    std::vector<Rect> touched;
    const Rect& p = s.layout.plot;
    const bool in_plot = s.ui.inside && p.contains(s.ui.mx, s.ui.my);

    // Readout: two lines at the bottom-left, SigPlot style.
    const Rect ro = s.layout.readout;
    fb.fill_rect(ro, th.bg);
    touched.push_back(ro);
    std::string line1 = "L=" + std::to_string(s.views.level()) + " " + cmode_short(s.set.cmode);
    std::string line2 = s.set.cmode == CMode::IR ? "(re/im)" : absc_label(s.set.absc);
    if (in_plot) {
        // SigPlot layout: "y: val dy: val L=n Ma" / "x: val dx: val (absc)", deltas from the marker.
        const double x = s.px_to_x(s.ui.mx), y = s.py_to_y(s.ui.my);
        std::string ys = "y: " + y_text(s, y), xs = "x: " + x_text(s, x, false);
        if (const auto z = s.content->z_at(x, y, s.set)) ys += "  z: " + format_g(*z, 9);
        if (s.ui.marker) {
            ys += "  dy: " + y_text(s, y - s.ui.marker->second);
            xs += "  dx: " + x_text(s, x - s.ui.marker->first, true);
        }
        line1 = ys + "  " + line1;
        line2 = xs + "  " + line2;
    }
    if (s.ui.prompt) {  // the prompt replaces the readout while it is open
        const auto kind = s.ui.prompt->kind;
        const char* axis = kind == PromptState::Kind::XRange ? "x" : kind == PromptState::Kind::YRange ? "y" : "z (colour)";
        line1 = std::string(axis) + " range (min max): " + s.ui.prompt->text + "_";
        line2 = s.ui.prompt->error.empty() ? "Enter to apply, Esc to cancel" : s.ui.prompt->error;
    }
    fb.text(kCW, ro.y + 4, line1, th.fg);
    fb.text(kCW, ro.y + 4 + kCH + 2, line2, th.fg);

    // Status, right-aligned on the first readout line.
    std::string status = s.message;
    const bool partial = s.result && !s.result->complete();
    if (rendering) {
        status = "rendering " + percent(progress) + " - click or Esc to stop";
    } else if (s.stopped && partial) {
        status = "stopped at " + percent(progress) + " - Space to resume";
    }
    if (!status.empty()) {
        fb.text(ro.right() - kCW - Framebuffer::text_width(status), ro.y + 4, status, th.fg);
    }
    static const std::string kHint = "?: help   middle-click: menu";
    fb.text(ro.right() - kCW - Framebuffer::text_width(kHint), ro.y + 4 + kCH + 2, kHint, th.dim);

    // Progress bar just above the data area.
    const Rect bar_area{p.x, p.y - 4, p.w, 2};
    if (rendering || (s.stopped && partial)) {
        const int w = static_cast<int>(std::lround(p.w * std::clamp(progress, 0.0, 1.0)));
        fb.fill_rect({bar_area.x, bar_area.y, w, bar_area.h}, rendering ? th.fg : th.dim);
        touched.push_back(bar_area);
    }

    if (s.ui.dragging) {
        const int x0 = std::clamp(std::min(s.ui.x0, s.ui.mx), p.x, p.right() - 1);
        const int x1 = std::clamp(std::max(s.ui.x0, s.ui.mx), p.x, p.right() - 1);
        const int y0 = std::clamp(std::min(s.ui.y0, s.ui.my), p.y, p.bottom() - 1);
        const int y1 = std::clamp(std::max(s.ui.y0, s.ui.my), p.y, p.bottom() - 1);
        fb.hline(x0, x1, y0, th.fg);
        fb.hline(x0, x1, y1, th.fg);
        fb.vline(x0, y0, y1, th.fg);
        fb.vline(x1, y0, y1, th.fg);
        touched.push_back({x0, y0, x1 - x0 + 1, 1});
        touched.push_back({x0, y1, x1 - x0 + 1, 1});
        touched.push_back({x0, y0, 1, y1 - y0 + 1});
        touched.push_back({x1, y0, 1, y1 - y0 + 1});
    } else if (in_plot && s.set.cross && !s.ui.menu && !s.ui.help) {
        fb.vline(s.ui.mx, p.y, p.bottom() - 1, th.fg, kCrossDot);
        fb.hline(p.x, p.right() - 1, s.ui.my, th.fg, kCrossDot);
        touched.push_back({s.ui.mx, p.y, 1, p.h});
        touched.push_back({p.x, s.ui.my, p.w, 1});
    }
    if (s.ui.marker) {  // small plus at the marker, with its value
        const Range X = s.xshown, Y = s.yshown;
        const double mxv = s.ui.marker->first, myv = s.ui.marker->second;
        const int px = p.x + static_cast<int>(std::floor((mxv - X.lo) / (X.hi - X.lo) * p.w));
        const double from_top = s.ydown() ? myv - Y.lo : Y.hi - myv;
        const int py = p.y + static_cast<int>(std::floor(from_top / (Y.hi - Y.lo) * p.h));
        if (p.contains(px, py)) {
            fb.hline(std::max(p.x, px - 4), std::min(p.right() - 1, px + 4), py, th.fg);
            fb.vline(px, std::max(p.y, py - 4), std::min(p.bottom() - 1, py + 4), th.fg);
            const std::string label = "x:" + format_g(mxv, 6) + " y:" + format_g(myv, 6);
            const int tw = Framebuffer::text_width(label);
            const int tx = std::min(px + 6, p.right() - tw), ty = std::max(p.y, py - kCH - 4);
            fb.fill_rect({tx, ty, tw, kCH}, th.bg);
            fb.text(tx, ty, label, th.fg);
            touched.push_back({px - 4, py - 4, 9, 9});
            touched.push_back({tx, ty, tw, kCH});
        }
    }
    for (const Rect& r : draw_menu(fb, s, th)) touched.push_back(r);
    if (s.ui.help) touched.push_back(draw_help(fb, th));
    return touched;
}

}  // namespace ssp
