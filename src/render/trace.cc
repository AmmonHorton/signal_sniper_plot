#include "render/trace.h"

#include <algorithm>
#include <cmath>

#include "core/scan.h"

namespace ssp {

Span bins_extent(const TraceBins& bins, int done) {
    Span s;
    done = std::min<int>(done, static_cast<int>(bins.cols.size()));
    for (int c = 0; c < done; ++c) {
        if (bins.cols[c].count) {
            s.add(bins.cols[c].lo);
            s.add(bins.cols[c].hi);
        }
    }
    return s;
}

TraceReducer::TraceReducer(Lod& lod, double xstart, double xdelta, const XView& view,
                           TraceBins& out)
    : lod_(lod), xstart_(xstart), xdelta_(xdelta), view_(view), out_(out),
      comp_(comp_of(view.cmode)) {
    prev_ = boundary(0);
    if (prev_ > 0) out_.before = {true, xstart_ + (prev_ - 1) * xdelta_, value(prev_ - 1)};
}

std::size_t TraceReducer::boundary(int c) const {
    // First sample index whose x >= the left edge of column c (c == width is the right edge).
    const std::size_t n = lod_.signal().n;
    const double col_dx = (view_.x1 - view_.x0) / view_.width;
    const double t = std::ceil((view_.x0 + c * col_dx - xstart_) / xdelta_);
    if (!(t > 0.0)) return 0;
    if (t >= static_cast<double>(n)) return n;
    return static_cast<std::size_t>(t);
}

double TraceReducer::value(std::size_t i) const {
    return display_value(view_.cmode, view_.units, sample_value(lod_.signal(), comp_, i));
}

void TraceReducer::run_to(int col_end, const CancelToken& ct) {
    const int w = static_cast<int>(out_.cols.size());
    col_end = std::min(col_end, w);
    for (; next_ < col_end; ++next_) {
        const std::size_t next = boundary(next_ + 1);
        if (next > prev_) {
            const Span s = lod_.query(comp_, prev_, next, ct);
            if (!s.empty()) {
                Bin& b = out_.cols[next_];
                b.lo = display_value(view_.cmode, view_.units, s.lo);
                b.hi = display_value(view_.cmode, view_.units, s.hi);
                b.first = value(prev_);
                b.last = value(next - 1);
                b.count = next - prev_;
                if (!out_.occ.empty() && b.count > 1) mark_rows(next_, prev_, next, b, ct);
            }
        }
        prev_ = next;
    }
    if (next_ == w && prev_ < lod_.signal().n) {
        out_.after = {true, xstart_ + prev_ * xdelta_, value(prev_)};
    }
}

Span TraceReducer::view_extent(const CancelToken& ct) const {
    const Span s = lod_.query(comp_, boundary(0), boundary(view_.width), ct);
    Span out;
    if (!s.empty()) {
        out.add(display_value(view_.cmode, view_.units, s.lo));
        out.add(display_value(view_.cmode, view_.units, s.hi));
    }
    return out;
}

void TraceReducer::mark_rows(int c, std::size_t i0, std::size_t i1, const Bin& b,
                             const CancelToken& ct) {
    const double y0 = out_.occ_y.lo, y1 = out_.occ_y.hi;
    const int h = view_.height;
    uint64_t* col = &out_.occ[static_cast<std::size_t>(c) * out_.occ_words];
    const int top = row_of(b.hi, y0, y1, h), bottom = row_of(b.lo, y0, y1, h);
    if (bottom < 0 || top >= h) return;  // whole column off-screen
    if (top == bottom) {                  // every sample lands in one row: no need to look
        col[top / 64] |= uint64_t{1} << (top % 64);
        return;
    }
    const double scale = h / (y1 - y0);
    for_each_value(lod_.signal(), comp_, i0, i1, ct, [&](double v) {
        v = display_value(view_.cmode, view_.units, v);
        const double r = std::floor((y1 - v) * scale);
        if (r >= 0.0 && r < h) {
            const int ri = static_cast<int>(r);
            col[ri / 64] |= uint64_t{1} << (ri % 64);
        }
    });
}

void reduce_trace(Lod& lod, double xstart, double xdelta, const XView& view, TraceBins& out,
                  const CancelToken& ct) {
    out = TraceBins{};
    out.cols.assign(std::max(1, view.width), Bin{});
    TraceReducer(lod, xstart, xdelta, view, out).run_to(view.width, ct);
}

void paint_trace(Framebuffer& fb, const TraceBins& bins, int done, const XView& view, Rect plot,
                 double y0, double y1, const TraceStyle& style) {
    // `before` is final once any column is published and `after` once all are; another
    // thread may still be writing either before that, so check `done` first.
    done = std::min<int>(done, static_cast<int>(bins.cols.size()));
    if (done <= 0 || plot.empty() || !(y1 > y0) || !(view.x1 > view.x0)) return;
    const bool all_done = done == static_cast<int>(bins.cols.size());
    const int thick = std::max(1, style.thickness);
    const bool lines = style.style != Style::Dots;
    const bool dots = style.style != Style::Lines;
    const int radius = thick;
    const bool exact_rows = !bins.occ.empty() && bins.occ_y.lo == y0 && bins.occ_y.hi == y1 &&
                            bins.occ_words * 64 >= plot.h;

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
            } else if (!lines && exact_rows) {
                const int x = plot.x + c - (thick - 1) / 2;
                for (int r = 0; r < plot.h; ++r) {
                    if (bins.has_row(c, r)) {
                        fb.fill_rect(Rect{x, plot.y + r - (thick - 1) / 2, thick, thick}.intersect(plot),
                                     style.color);
                    }
                }
            } else if (!lines) {
                span(c, b.lo, b.hi);  // no row data (e.g. headless callers that skip it)
            }
        }
    }
    if (lines && have_prev && all_done && bins.after.valid) {
        fb.line(prev_x, prev_y, px_of_x(bins.after.x), py(bins.after.v), style.color, thick, plot);
    }
}

}  // namespace ssp
