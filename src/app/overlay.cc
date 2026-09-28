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
    "  any click        stop a render in progress",
    "",
    "Keys",
    "  1 2 3 4 6 7      Magnitude Phase Real Imag 10log 20log",
    "  Home             unzoom all",
    "  Space / Esc      resume / stop rendering",
    "  i                index x-axis",
    "  g  l             grid, legend",
    "  m                menu (also: x/y range, phase units, traces)",
    "  Ctrl-S           save PNG",
    "  q                quit",
    "  ?                this help (any key closes)",
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
    std::string line2 = s.set.index ? "(indx)" : "(absc)";
    if (in_plot) {
        const double x = s.px_to_x(s.ui.mx);
        const std::string xs = s.set.index ? std::to_string(std::llround(x)) : format_g(x, 9);
        line1 = "y: " + format_g(s.py_to_y(s.ui.my), 9) + "  " + line1;
        line2 = "x: " + xs + "  " + line2;
    }
    if (s.ui.prompt) {  // the prompt replaces the readout while it is open
        const bool x = s.ui.prompt->kind == PromptState::Kind::XRange;
        line1 = std::string(x ? "x" : "y") + " range (min max): " + s.ui.prompt->text + "_";
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
    for (const Rect& r : draw_menu(fb, s, th)) touched.push_back(r);
    if (s.ui.help) touched.push_back(draw_help(fb, th));
    return touched;
}

}  // namespace ssp
