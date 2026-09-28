#include "app/raster_content.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "core/dispatch.h"
#include "core/lod.h"
#include "core/scan.h"
#include "render/colormap.h"

namespace ssp {
namespace {

constexpr auto kPublishEvery = std::chrono::milliseconds(8);

/// Which cells each of `n` pixels covers along one axis showing [lo, hi): the cells whose
/// centres fall inside the pixel, or (zoomed in past one cell per pixel) the cell under the
/// pixel's centre.
struct PixelCells {
    std::vector<std::size_t> edge;   ///< n + 1 cell boundaries: pixel i has [edge[i], edge[i+1])
    std::vector<long long> nearest;  ///< Cell under pixel i's centre, -1 if none.

    /// Cells [first, last) of pixel i (possibly empty).
    std::pair<std::size_t, std::size_t> of(int i) const {
        if (edge[i + 1] > edge[i]) return {edge[i], edge[i + 1]};
        if (nearest[i] >= 0) return {std::size_t(nearest[i]), std::size_t(nearest[i]) + 1};
        return {0, 0};
    }
};

PixelCells pixel_cells(double lo, double hi, int n, Range axis, std::size_t count) {
    PixelCells p;
    p.edge.resize(n + 1);
    p.nearest.resize(n);
    const double d = (hi - lo) / n;
    for (int i = 0; i <= n; ++i) {
        const double t = std::ceil((lo + i * d - axis.lo) / axis.hi - 0.5);
        p.edge[i] = !(t > 0.0) ? 0 : t >= double(count) ? count : std::size_t(t);
    }
    for (int i = 0; i < n; ++i) {
        const double c = std::floor((lo + (i + 0.5) * d - axis.lo) / axis.hi);
        p.nearest[i] = (c >= 0.0 && c < double(count)) ? static_cast<long long>(c) : -1;
    }
    return p;
}

/// Running reduction of one pixel.
struct Acc {
    double v = 0.0;
    uint32_t n = 0;

    void add(Reduce mode, double x) {
        if (!(x == x)) return;  // NaN samples don't count
        switch (mode) {
            case Reduce::Max:    v = n ? std::max(v, x) : x; break;
            case Reduce::Min:    v = n ? std::min(v, x) : x; break;
            case Reduce::MaxAbs: v = n ? std::max(v, std::abs(x)) : std::abs(x); break;
            case Reduce::First:  if (!n) v = x; break;
            case Reduce::Mean:   v += x; break;
        }
        ++n;
    }
    double result(Reduce mode) const { return mode == Reduce::Mean ? v / n : v; }
};

void check_range(const std::optional<Range>& r, const char* what) {
    if (r && !(std::isfinite(r->lo) && std::isfinite(r->hi) && r->hi > r->lo)) {
        throw std::invalid_argument(std::string(what) + " must be finite with hi > lo");
    }
}

}  // namespace

RasterContent::RasterContent(Signal data, RasterOptions o) : sig_(std::move(data)), opts_(std::move(o)) {
    if (opts_.subsize == 0) throw std::invalid_argument("raster: subsize (samples per frame) must be > 0");
    if (sig_.n == 0) throw std::invalid_argument("raster: no data");
    if (!sig_.data) throw std::invalid_argument("raster: data is null");
    if (sig_.stride < 1) throw std::invalid_argument("raster: stride must be >= 1");
    if (!(std::isfinite(sig_.xdelta) && sig_.xdelta > 0.0) || !std::isfinite(sig_.xstart)) {
        throw std::invalid_argument("raster: xdelta must be finite and > 0");
    }
    if (!(std::isfinite(opts_.ydelta) && opts_.ydelta > 0.0) || !std::isfinite(opts_.ystart)) {
        throw std::invalid_argument("raster: ydelta must be finite and > 0");
    }
    if (opts_.frame_stride < 0) throw std::invalid_argument("raster: frame_stride must be >= 0");
    if (opts_.frame_stride == 0) opts_.frame_stride = static_cast<std::ptrdiff_t>(opts_.subsize) * sig_.stride;
    check_range(opts_.xrange, "xrange");
    check_range(opts_.yrange, "yrange");
    check_range(opts_.zrange, "zrange");
    rows_ = (sig_.n + opts_.subsize - 1) / opts_.subsize;
    initial_cmode_ = opts_.cmode == CMode::Auto ? (sig_.complex ? CMode::Mag : CMode::Real) : opts_.cmode;
    if (initial_cmode_ == CMode::IR) throw std::invalid_argument("raster: IR mode is not available");
}

Settings RasterContent::initial_settings() const {
    Settings s;
    s.title = opts_.title;
    s.cmode = initial_cmode_;
    s.phunits = opts_.phunits;
    s.index = opts_.index;
    s.grid = opts_.grid;
    s.legend = false;
    s.cmap = opts_.cmap;
    s.reduce = opts_.reduce;
    s.zfixed = opts_.zrange;
    return s;
}

View RasterContent::home(const Settings& s) const {
    const bool original_axes = s.index == opts_.index;
    const Range cx = col_axis(s.index), ry = row_axis(s.index);
    View v;
    v.x = (original_axes && opts_.xrange) ? *opts_.xrange : Range{cx.lo, cx.lo + cols() * cx.hi};
    v.y = (original_axes && opts_.yrange) ? *opts_.yrange : Range{ry.lo, ry.lo + rows_ * ry.hi};
    return v;
}

std::size_t RasterContent::cols_in_row(std::size_t r) const {
    const std::size_t first = r * opts_.subsize;
    return first >= sig_.n ? 0 : std::min(opts_.subsize, sig_.n - first);
}

Signal RasterContent::row_signal(std::size_t r, std::size_t c0, std::size_t c1) const {
    const std::size_t elem = dtype_size(sig_.dtype) * (sig_.complex ? 2 : 1);
    const std::ptrdiff_t off = static_cast<std::ptrdiff_t>(r) * opts_.frame_stride +
                               static_cast<std::ptrdiff_t>(c0) * sig_.stride;
    Signal s = sig_;
    s.data = static_cast<const char*>(sig_.data) + off * static_cast<std::ptrdiff_t>(elem);
    s.n = c1 > c0 ? c1 - c0 : 0;
    s.xstart = sig_.xstart + c0 * sig_.xdelta;
    return s;
}

Signal RasterContent::col_signal(std::size_t c, std::size_t r0, std::size_t r1) const {
    Signal s = row_signal(r0, c, c + 1);
    s.n = r1 > r0 ? r1 - r0 : 0;
    s.stride = opts_.frame_stride;
    s.xstart = opts_.ystart + r0 * opts_.ydelta;
    s.xdelta = opts_.ydelta;
    return s;
}

std::optional<std::size_t> RasterContent::col_at(double x, bool index) const {
    const Range a = col_axis(index);
    const double c = std::floor((x - a.lo) / a.hi);
    if (!(c >= 0.0 && c < double(cols()))) return std::nullopt;
    return static_cast<std::size_t>(c);
}

std::optional<std::size_t> RasterContent::row_at(double y, bool index) const {
    const Range a = row_axis(index);
    const double r = std::floor((y - a.lo) / a.hi);
    if (!(r >= 0.0 && r < double(rows_))) return std::nullopt;
    return static_cast<std::size_t>(r);
}

namespace {
SampleRange overlapping(Range v, Range axis, std::size_t count) {
    auto clamp_cell = [&](double t) -> std::size_t {
        return !(t > 0.0) ? 0 : t >= double(count) ? count : std::size_t(t);
    };
    const std::size_t first = clamp_cell(std::floor((v.lo - axis.lo) / axis.hi));
    const std::size_t last = clamp_cell(std::ceil((v.hi - axis.lo) / axis.hi));
    return {first, std::max(first, last)};
}
}  // namespace

SampleRange RasterContent::cols_in(Range x, bool index) const { return overlapping(x, col_axis(index), cols()); }
SampleRange RasterContent::rows_in(Range y, bool index) const { return overlapping(y, row_axis(index), rows_); }

std::shared_ptr<ReduceResult> RasterContent::new_result(const ReduceRequest& req) const {
    auto r = std::make_shared<ReduceResult>();
    r->req = req;
    r->view = req.view;
    r->view.width = std::max(1, req.view.width);
    r->view.height = std::max(1, req.view.height);
    r->total = r->view.height;
    r->x = Range{req.view.x0, req.view.x1};
    r->y = req.y ? *req.y : home(initial_settings()).y;
    r->zimg.assign(static_cast<std::size_t>(r->view.width) * r->view.height,
                   std::numeric_limits<float>::quiet_NaN());
    return r;
}

bool RasterContent::reusable(const ReduceResult& r, const ReduceRequest& req) const {
    if (!r.complete()) return false;
    const ReduceRequest& q = r.req;
    const XView &a = q.view, &b = req.view;
    const bool same_y = q.y.has_value() == req.y.has_value() &&
                        (!q.y || (q.y->lo == req.y->lo && q.y->hi == req.y->hi));
    return a.x0 == b.x0 && a.x1 == b.x1 && a.width == b.width && a.height == b.height &&
           a.cmode == b.cmode && a.units == b.units && q.index == req.index &&
           q.reduce == req.reduce && same_y;
}

void RasterContent::reduce(ReduceResult& r, const CancelToken& ct,
                           const std::function<void()>& progress) const {
    const int w = r.view.width, h = r.view.height;
    const Range x = *r.x, y = *r.y;
    const PixelCells px_cols = pixel_cells(x.lo, x.hi, w, col_axis(r.req.index), cols());
    // Frame 0 is at the top, so pixel rows run from y.lo downwards.
    const PixelCells px_rows = pixel_cells(y.lo, y.hi, h, row_axis(r.req.index), rows_);
    const Comp comp = comp_of(r.view.cmode);
    const Reduce mode = r.req.reduce;
    // Max/Min/First commute with the (monotonic) display transform, so it is applied once per
    // pixel; Mean and MaxAbs need it per sample.
    const bool per_sample = mode == Reduce::Mean || mode == Reduce::MaxAbs;
    auto disp = [&](double v) { return display_value(r.view.cmode, r.view.units, v); };

    std::size_t c_lo = cols(), c_hi = 0;  // columns anywhere in view
    for (int i = 0; i < w; ++i) {
        const auto [a, b] = px_cols.of(i);
        if (a < b) {
            c_lo = std::min(c_lo, a);
            c_hi = std::max(c_hi, b);
        }
    }

    std::vector<double> vals;
    std::vector<Acc> acc(w);
    std::pair<std::size_t, std::size_t> prev_rows{1, 0};  // impossible: forces the first compute
    auto last_publish = std::chrono::steady_clock::now();

    for (int py = 0; py < h; ++py) {
        float* out = &r.zimg[static_cast<std::size_t>(py) * w];
        const auto rows = px_rows.of(py);
        if (rows == prev_rows && py > 0) {  // zoomed in: same frame rows as the pixel row above
            std::copy(out - w, out, out);
        } else if (rows.first < rows.second && c_lo < c_hi) {
            std::fill(acc.begin(), acc.end(), Acc{});
            for (std::size_t fr = rows.first; fr < rows.second; ++fr) {
                const std::size_t ce = std::min(c_hi, cols_in_row(fr));
                if (c_lo >= ce) continue;
                vals.resize(ce - c_lo);
                std::size_t k = 0;
                for_each_value(row_signal(fr, c_lo, ce), comp, 0, ce - c_lo, ct,
                               [&](double v) { vals[k++] = per_sample ? disp(v) : v; });
                for (int i = 0; i < w; ++i) {
                    auto [a, b] = px_cols.of(i);
                    b = std::min(b, ce);
                    for (std::size_t c = a; c < b; ++c) acc[i].add(mode, vals[c - c_lo]);
                }
            }
            for (int i = 0; i < w; ++i) {
                if (acc[i].n) {
                    const double v = acc[i].result(mode);
                    out[i] = static_cast<float>(per_sample ? v : disp(v));
                }
            }
        }
        prev_rows = rows;
        const auto now = std::chrono::steady_clock::now();
        if (py + 1 == h || now - last_publish >= kPublishEvery) {
            r.done.store(py + 1, std::memory_order_release);
            last_publish = now;
            if (progress) progress();
        }
        ct.check();
    }
}

Span RasterContent::extent(const ReduceResult& r, int done) const {
    Span s;
    const std::size_t n = static_cast<std::size_t>(std::min(done, r.view.height)) * r.view.width;
    for (std::size_t k = 0; k < n; ++k) s.add(r.zimg[k]);
    return s;
}

void RasterContent::paint(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const {
    const int w = r.view.width, h = r.view.height;
    if (a.done <= 0 || w != a.plot.w || h != a.plot.h || !(a.z.hi > a.z.lo)) return;
    const ColorLut& lut = colormap_lut(a.set.cmap);
    const int rows = std::min(a.done, h);
    for (int py = 0; py < rows; ++py) {
        uint32_t* dst = fb.row(a.plot.y + py) + a.plot.x;
        const float* z = &r.zimg[static_cast<std::size_t>(py) * w];
        for (int i = 0; i < w; ++i) {
            const int k = lut_index(z[i], a.z.lo, a.z.hi);
            if (k >= 0) dst[i] = lut[k];
        }
    }
}

std::optional<double> RasterContent::z_at(double x, double y, const Settings& s) const {
    const auto c = col_at(x, s.index), r = row_at(y, s.index);
    if (!c || !r || *c >= cols_in_row(*r)) return std::nullopt;
    const double v = sample_value(row_signal(*r, *c, *c + 1), comp_of(s.cmode), 0);
    return display_value(s.cmode, s.phunits, v);
}

}  // namespace ssp
