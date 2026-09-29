#include "ssp/render/frame.h"

#include <algorithm>
#include <cmath>

namespace ssp {
namespace {

constexpr int kCW = Framebuffer::kCharW;
constexpr int kCH = Framebuffer::kCharH;
constexpr int kTop = 2 * kCH + 4;          // title row (+ room for a progress bar)
constexpr int kXLabels = kCH + 6;          // x tick labels
constexpr int kReadout = 2 * kCH + 6;      // two readout lines
constexpr int kRight = 16;
constexpr int kBarGap = 10, kBarW = 12;                                       // colour bar
constexpr int kBarRight = kBarGap + kBarW + 6 + kMaxTickLabelChars * kCW + 6;  // + its labels
constexpr int kTickLen = 4;
constexpr int kGridDot = 3;

}  // namespace

Layout compute_layout(int window_w, int window_h, bool colorbar) {
    constexpr int kLeft = kMaxTickLabelChars * kCW + 12;
    const int right = colorbar ? kRight + kBarRight : kRight;
    Layout l;
    l.plot = {kLeft, kTop, std::max(1, window_w - kLeft - right),
              std::max(1, window_h - kTop - kXLabels - kReadout)};
    l.title = {0, 0, window_w, kTop};
    l.xlabels = {0, l.plot.bottom(), window_w, kXLabels};
    l.readout = {0, l.xlabels.bottom(), window_w, kReadout};
    return l;
}

int xdivisions(int plot_w) { return std::clamp(plot_w / 110, 2, 10); }
int ydivisions(int plot_h) { return std::clamp(plot_h / 60, 2, 10); }

void draw_axes(Framebuffer& fb, const Layout& l, double x0, double x1, const AxisTicks& xt,
               double y0, double y1, const AxisTicks& yt, bool grid, const Theme& th,
               bool ydown) {
    const Rect& p = l.plot;

    // x: grid, tick marks, labels (skipping any that would overlap the previous one or the note).
    const int note_x = xt.note.empty() ? p.right() + kRight
                                       : p.right() - Framebuffer::text_width(xt.note);
    int last_label_end = -1;
    for (std::size_t i = 0; i < xt.values.size(); ++i) {
        // Ticks are inside [x0, x1]; the one at x1 maps to the last column, not past it.
        const int x = std::min(p.right() - 1,
                               p.x + static_cast<int>(std::floor((xt.values[i] - x0) / (x1 - x0) * p.w)));
        if (x < p.x) continue;
        if (grid) fb.vline(x, p.y, p.bottom() - 1, th.grid, kGridDot);
        fb.vline(x, p.bottom() - kTickLen, p.bottom() - 1, th.fg);
        fb.vline(x, p.y, p.y + kTickLen - 1, th.fg);
        const std::string& s = xt.labels[i];
        const int tw = Framebuffer::text_width(s);
        const int tx = std::clamp(x - tw / 2, 0, fb.width() - tw);
        if (tx > last_label_end + kCW && tx + tw < note_x - kCW) {
            fb.text(tx, l.xlabels.y + 4, s, th.fg);
            last_label_end = tx + tw;
        }
    }
    if (!xt.note.empty()) fb.text(note_x, l.xlabels.y + 4, xt.note, th.fg);

    // y
    for (std::size_t i = 0; i < yt.values.size(); ++i) {
        const double from_top = ydown ? yt.values[i] - y0 : y1 - yt.values[i];
        const int y = std::min(p.bottom() - 1, p.y + static_cast<int>(std::floor(from_top / (y1 - y0) * p.h)));
        if (y < p.y) continue;
        if (grid) fb.hline(p.x, p.right() - 1, y, th.grid, kGridDot);
        fb.hline(p.x, p.x + kTickLen - 1, y, th.fg);
        fb.hline(p.right() - kTickLen, p.right() - 1, y, th.fg);
        const std::string& s = yt.labels[i];
        const int ty = std::clamp(y - kCH / 2, p.y - kCH / 2, p.bottom() - kCH / 2 - 1);
        fb.text(p.x - 6 - Framebuffer::text_width(s), ty, s, th.fg);
    }
    if (!yt.note.empty()) fb.text(2, p.y - kCH - 2, yt.note, th.fg);

    fb.rect_outline({p.x - 1, p.y - 1, p.w + 2, p.h + 2}, th.fg);
}

void draw_colorbar(Framebuffer& fb, const Layout& l, double z0, double z1,
                   const std::array<uint32_t, 256>& lut, const Theme& th) {
    const Rect& p = l.plot;
    const Rect bar{p.right() + kBarGap, p.y, kBarW, p.h};
    for (int y = 0; y < bar.h; ++y) {
        const int k = std::clamp(static_cast<int>((1.0 - (y + 0.5) / bar.h) * 256.0), 0, 255);
        fb.hline(bar.x, bar.right() - 1, bar.y + y, lut[k]);
    }
    fb.rect_outline({bar.x - 1, bar.y - 1, bar.w + 2, bar.h + 2}, th.fg);
    if (!(z1 > z0)) return;
    const AxisTicks zt = make_ticks(z0, z1, ydivisions(bar.h));
    for (std::size_t i = 0; i < zt.values.size(); ++i) {
        const int y = std::min(bar.bottom() - 1,
                               bar.y + static_cast<int>(std::floor((z1 - zt.values[i]) / (z1 - z0) * bar.h)));
        fb.hline(bar.right() + 1, bar.right() + 3, y, th.fg);
        fb.text(bar.right() + 6, std::clamp(y - kCH / 2, bar.y - kCH / 2, bar.bottom() - kCH / 2), zt.labels[i], th.fg);
    }
    if (!zt.note.empty()) fb.text(bar.x, bar.bottom() + 4, zt.note, th.fg);
}

void draw_title(Framebuffer& fb, const Layout& l, const std::string& title, const Theme& th) {
    if (title.empty()) return;
    const int tw = Framebuffer::text_width(title);
    fb.text(std::max(0, l.title.x + (l.title.w - tw) / 2), l.title.y + 6, title, th.fg);
}

std::vector<Rect> draw_legend(Framebuffer& fb, const Rect& plot,
                              const std::vector<LegendEntry>& entries, const Theme& th) {
    std::vector<Rect> hits;
    if (entries.empty()) return hits;
    constexpr int kSwatch = 18, kRow = kCH + 2, kPad = 4;
    std::size_t chars = 0;
    for (const auto& e : entries) chars = std::max(chars, e.name.size());
    const int w = kPad + kSwatch + kCW + static_cast<int>(chars) * kCW + kPad;
    const int h = kPad + static_cast<int>(entries.size()) * kRow + kPad;
    const Rect box{plot.right() - w - 6, plot.y + 6, w, h};

    fb.fill_rect(box, th.bg);
    fb.rect_outline(box, th.fg);
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const auto& e = entries[i];
        const int y = box.y + kPad + static_cast<int>(i) * kRow;
        const uint32_t c = e.visible ? e.color : th.dim;
        const int mid = y + kRow / 2;
        if (e.style != Style::Dots) fb.hline(box.x + kPad, box.x + kPad + kSwatch - 1, mid, c);
        if (e.style != Style::Lines) fb.disk(box.x + kPad + kSwatch / 2, mid, 2, c, box);
        fb.text(box.x + kPad + kSwatch + kCW, y + 1, e.name, e.visible ? th.fg : th.dim);
        hits.push_back({box.x, y, box.w, kRow});
    }
    return hits;
}

const char* cmode_name(CMode m) {
    switch (m) {
        case CMode::Auto:  return "Auto";
        case CMode::Mag:   return "Magnitude";
        case CMode::Phase: return "Phase";
        case CMode::Real:  return "Real";
        case CMode::Imag:  return "Imaginary";
        case CMode::IR:    return "Imag vs Real";
        case CMode::Log10: return "10*log10";
        case CMode::Log20: return "20*log10";
    }
    return "?";
}

}  // namespace ssp
