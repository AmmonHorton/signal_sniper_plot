#include "render/trace.h"

#include <algorithm>
#include <cmath>

namespace ssp {

void reduce_trace(Lod& lod, double xstart, double xdelta, const XView& view, TraceBins& out,
                  const CancelToken& ct, std::atomic<int>* done) {
    const std::size_t n = lod.signal().n;
    const int w = std::max(1, view.width);
    const Comp comp = comp_of(view.cmode);
    const double col_dx = (view.x1 - view.x0) / w;

    // First sample index whose x >= the left edge of column c (c == w is the right edge).
    auto boundary = [&](int c) -> std::size_t {
        const double t = std::ceil((view.x0 + c * col_dx - xstart) / xdelta);
        if (!(t > 0.0)) return 0;
        if (t >= static_cast<double>(n)) return n;
        return static_cast<std::size_t>(t);
    };
    auto value = [&](std::size_t i) {
        return display_value(view.cmode, view.units, sample_value(lod.signal(), comp, i));
    };

    out.cols.assign(w, Bin{});
    out.extent = {};
    out.before = out.after = {};

    std::size_t prev = boundary(0);
    if (prev > 0) out.before = {true, xstart + (prev - 1) * xdelta, value(prev - 1)};

    for (int c = 0; c < w; ++c) {
        const std::size_t next = boundary(c + 1);
        if (next > prev) {
            const Span s = lod.query(comp, prev, next, ct);
            if (!s.empty()) {
                Bin& b = out.cols[c];
                b.lo = display_value(view.cmode, view.units, s.lo);
                b.hi = display_value(view.cmode, view.units, s.hi);
                b.first = value(prev);
                b.last = value(next - 1);
                b.count = next - prev;
                out.extent.add(b.lo);
                out.extent.add(b.hi);
            }
        }
        prev = next;
        if (done) done->store(c + 1, std::memory_order_release);
    }
    if (prev < n) out.after = {true, xstart + prev * xdelta, value(prev)};
}

void paint_trace(Framebuffer& fb, const TraceBins& bins, int done, const XView& view, Rect plot,
                 double y0, double y1, const TraceStyle& style) {
    if (plot.empty() || !(y1 > y0) || !(view.x1 > view.x0)) return;
    done = std::min<int>(done, static_cast<int>(bins.cols.size()));
    const int thick = std::max(1, style.thickness);
    const bool lines = style.style != Style::Dots;
    const bool dots = style.style != Style::Lines;
    const int radius = thick;

    auto py = [&](double v) { return plot.y + (y1 - v) / (y1 - y0) * plot.h; };
    auto px_of_x = [&](double x) {
        return plot.x + (x - view.x0) / (view.x1 - view.x0) * bins.cols.size();
    };

    // Vertical run covering [lo, hi] in column c, if any of it is visible.
    auto span = [&](int c, double lo, double hi) {
        const double top = py(hi), bot = py(lo);
        if (bot < plot.y || top >= plot.bottom()) return;
        const int r0 = std::max(plot.y, static_cast<int>(std::floor(top)));
        const int r1 = std::min(plot.bottom() - 1, static_cast<int>(std::floor(bot)));
        const int x = plot.x + c - (thick - 1) / 2;
        fb.fill_rect(Rect{x, r0, thick, r1 - r0 + 1}.intersect(plot), style.color);
    };

    bool have_prev = false;
    double prev_x = 0.0, prev_y = 0.0;
    if (bins.before.valid) {
        have_prev = true;
        prev_x = px_of_x(bins.before.x);
        prev_y = py(bins.before.v);
    }

    for (int c = 0; c < done; ++c) {
        const Bin& b = bins.cols[c];
        if (b.count == 0) continue;
        const double cx = plot.x + c + 0.5;
        if (lines) {
            if (have_prev) fb.line(prev_x, prev_y, cx, py(b.first), style.color, thick, plot);
            span(c, b.lo, b.hi);
            have_prev = true;
            prev_x = cx;
            prev_y = py(b.last);
        }
        if (dots) {
            if (b.count == 1) {
                const double y = py(b.first);
                if (y >= plot.y && y < plot.bottom()) {
                    fb.disk(plot.x + c, static_cast<int>(std::floor(y)), radius, style.color, plot);
                }
            } else if (!lines) {
                span(c, b.lo, b.hi);
            }
        }
    }
    if (lines && have_prev && bins.after.valid && done == static_cast<int>(bins.cols.size())) {
        fb.line(prev_x, prev_y, px_of_x(bins.after.x), py(bins.after.v), style.color, thick, plot);
    }
}

}  // namespace ssp
