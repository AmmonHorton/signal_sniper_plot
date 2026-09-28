#include "render/framebuffer.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace ssp {
namespace {

constexpr uint8_t kFont[][Framebuffer::kCharH] = {
#include "render/font_6x13.inc"
};
constexpr int kFirstGlyph = 32;
constexpr int kLastGlyph = 126;

/// Liang–Barsky: clip segment p0→p1 to [xmin,xmax]×[ymin,ymax]. False if fully outside.
bool clip_segment(double& x0, double& y0, double& x1, double& y1,
                  double xmin, double xmax, double ymin, double ymax) {
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(x1) || !std::isfinite(y1)) {
        return false;
    }
    const double dx = x1 - x0, dy = y1 - y0;
    double t0 = 0.0, t1 = 1.0;
    const double p[4] = {-dx, dx, -dy, dy};
    const double q[4] = {x0 - xmin, xmax - x0, y0 - ymin, ymax - y0};
    for (int i = 0; i < 4; ++i) {
        if (p[i] == 0.0) {
            if (q[i] < 0.0) return false;
            continue;
        }
        const double t = q[i] / p[i];
        if (p[i] < 0.0) {
            if (t > t1) return false;
            t0 = std::max(t0, t);
        } else {
            if (t < t0) return false;
            t1 = std::min(t1, t);
        }
    }
    const double nx0 = x0 + t0 * dx, ny0 = y0 + t0 * dy;
    x1 = x0 + t1 * dx;
    y1 = y0 + t1 * dy;
    x0 = nx0;
    y0 = ny0;
    return true;
}

}  // namespace

Rect Rect::intersect(const Rect& o) const {
    const int nx = std::max(x, o.x), ny = std::max(y, o.y);
    const int nr = std::min(right(), o.right()), nb = std::min(bottom(), o.bottom());
    return {nx, ny, std::max(0, nr - nx), std::max(0, nb - ny)};
}

Framebuffer::Framebuffer(int w, int h, uint32_t fill)
    : w_(std::max(0, w)), h_(std::max(0, h)), px_(static_cast<std::size_t>(w_) * h_, fill) {}

void Framebuffer::fill_rect(Rect r, uint32_t c) {
    r = r.intersect(bounds());
    for (int y = r.y; y < r.bottom(); ++y) std::fill(row(y) + r.x, row(y) + r.right(), c);
}

void Framebuffer::rect_outline(Rect r, uint32_t c) {
    if (r.empty()) return;
    hline(r.x, r.right() - 1, r.y, c);
    hline(r.x, r.right() - 1, r.bottom() - 1, c);
    vline(r.x, r.y, r.bottom() - 1, c);
    vline(r.right() - 1, r.y, r.bottom() - 1, c);
}

void Framebuffer::hline(int x0, int x1, int y, uint32_t c, int dot) {
    if (y < 0 || y >= h_) return;
    if (x0 > x1) std::swap(x0, x1);
    x0 = std::max(x0, 0);
    x1 = std::min(x1, w_ - 1);
    for (int x = x0; x <= x1; ++x) {
        if (dot <= 0 || x % dot == 0) row(y)[x] = c;
    }
}

void Framebuffer::vline(int x, int y0, int y1, uint32_t c, int dot) {
    if (x < 0 || x >= w_) return;
    if (y0 > y1) std::swap(y0, y1);
    y0 = std::max(y0, 0);
    y1 = std::min(y1, h_ - 1);
    for (int y = y0; y <= y1; ++y) {
        if (dot <= 0 || y % dot == 0) row(y)[x] = c;
    }
}

void Framebuffer::brush(int x, int y, uint32_t c, int thick, const Rect& clip) {
    const int lo = -(thick - 1) / 2, hi = thick / 2;
    for (int dy = lo; dy <= hi; ++dy) {
        for (int dx = lo; dx <= hi; ++dx) {
            if (clip.contains(x + dx, y + dy)) row(y + dy)[x + dx] = c;
        }
    }
}

void Framebuffer::line(double x0, double y0, double x1, double y1, uint32_t c, int thick,
                       Rect clip) {
    clip = clip.intersect(bounds());
    if (clip.empty()) return;
    // Pixels cover [p, p+1); keep clipped endpoints strictly inside the last row/column.
    constexpr double kIn = 1e-9;
    if (!clip_segment(x0, y0, x1, y1, clip.x, clip.right() - kIn, clip.y, clip.bottom() - kIn)) {
        return;
    }
    int ix0 = static_cast<int>(std::floor(x0)), iy0 = static_cast<int>(std::floor(y0));
    const int ix1 = static_cast<int>(std::floor(x1)), iy1 = static_cast<int>(std::floor(y1));
    const int dx = std::abs(ix1 - ix0), sx = ix0 < ix1 ? 1 : -1;
    const int dy = -std::abs(iy1 - iy0), sy = iy0 < iy1 ? 1 : -1;
    int err = dx + dy;
    thick = std::max(1, thick);
    for (;;) {
        brush(ix0, iy0, c, thick, clip);
        if (ix0 == ix1 && iy0 == iy1) break;
        const int e2 = 2 * err;
        if (e2 >= dy) { err += dy; ix0 += sx; }
        if (e2 <= dx) { err += dx; iy0 += sy; }
    }
}

void Framebuffer::disk(int cx, int cy, int r, uint32_t c, Rect clip) {
    clip = clip.intersect(bounds());
    for (int dy = -r; dy <= r; ++dy) {
        for (int dx = -r; dx <= r; ++dx) {
            if (dx * dx + dy * dy <= r * r + r && clip.contains(cx + dx, cy + dy)) {
                row(cy + dy)[cx + dx] = c;
            }
        }
    }
}

void Framebuffer::text(int x, int y, std::string_view s, uint32_t c) {
    for (char ch : s) {
        int g = static_cast<unsigned char>(ch);
        if (g < kFirstGlyph || g > kLastGlyph) g = '?';
        const uint8_t* rows = kFont[g - kFirstGlyph];
        for (int ry = 0; ry < kCharH; ++ry) {
            const int py = y + ry;
            if (py < 0 || py >= h_ || rows[ry] == 0) continue;
            for (int rx = 0; rx < kCharW; ++rx) {
                const int px = x + rx;
                if (px >= 0 && px < w_ && (rows[ry] >> (kCharW - 1 - rx)) & 1) row(py)[px] = c;
            }
        }
        x += kCharW;
    }
}

}  // namespace ssp
