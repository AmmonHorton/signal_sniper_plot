#include "render/trace_plot.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ssp {
namespace {

constexpr double kYPad = 0.02;  // SigPlot pads autoscaled y by 2% each side

void check_range(const std::optional<Range>& r, const char* what) {
    if (r && !(std::isfinite(r->lo) && std::isfinite(r->hi) && r->hi > r->lo)) {
        throw std::invalid_argument(std::string(what) + " must be finite with hi > lo");
    }
}

Range widen_if_flat(Range r) {
    if (!(r.hi > r.lo)) return {r.lo - 1.0, r.hi + 1.0};
    return r;
}

}  // namespace

TracePlot::TracePlot(std::vector<Signal> signals, PlotOptions options)
    : sigs_(std::move(signals)), opts_(std::move(options)) {
    if (sigs_.empty()) throw std::invalid_argument("plot needs at least one signal");
    check_range(opts_.xrange, "xrange");
    check_range(opts_.yrange, "yrange");

    bool any_complex = false;
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        Signal& s = sigs_[i];
        const std::string who = "signal " + std::to_string(i) + ": ";
        if (s.n > 0 && !s.data) throw std::invalid_argument(who + "data is null");
        if (s.stride < 1) throw std::invalid_argument(who + "stride must be >= 1");
        if (!opts_.index && !(std::isfinite(s.xdelta) && s.xdelta > 0.0)) {
            throw std::invalid_argument(who + "xdelta must be finite and > 0");
        }
        if (!std::isfinite(s.xstart)) throw std::invalid_argument(who + "xstart must be finite");
        show_legend_ = show_legend_ || !s.name.empty() || sigs_.size() > 1;
        if (s.name.empty()) s.name = "Trace " + std::to_string(i);
        any_complex = any_complex || s.complex;
        lods_.push_back(std::make_unique<Lod>(s));
    }

    cmode_ = opts_.cmode;
    if (cmode_ == CMode::Auto) cmode_ = any_complex ? CMode::Mag : CMode::Real;
    if (cmode_ == CMode::IR) throw std::invalid_argument("IR (imag vs real) mode is not implemented yet");
    bins_.resize(sigs_.size());
}

Range TracePlot::autoscale_x() const {
    if (opts_.xrange) return *opts_.xrange;
    Span s;
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        if (!sigs_[i].visible || sigs_[i].n == 0) continue;
        s.add(xstart(i));
        s.add(xstart(i) + (sigs_[i].n - 1) * xdelta(i));
    }
    if (s.empty()) return {-1.0, 1.0};
    return widen_if_flat({s.lo, s.hi});
}

XView TracePlot::make_view(Range x, int width) const {
    return {x.lo, x.hi, std::max(1, width), cmode_, opts_.phunits};
}

void TracePlot::reduce(const XView& view, const CancelToken& ct) {
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        if (sigs_[i].visible) reduce_trace(*lods_[i], xstart(i), xdelta(i), view, bins_[i], ct);
    }
}

Range TracePlot::autoscale_y() const {
    if (opts_.yrange) return *opts_.yrange;
    Span s;
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        if (sigs_[i].visible) s.merge(bins_[i].extent);
    }
    if (s.empty()) return {-1.0, 1.0};
    const Range r = widen_if_flat({s.lo, s.hi});
    const double pad = kYPad * (r.hi - r.lo);
    return {r.lo - pad, r.hi + pad};
}

void TracePlot::paint(Framebuffer& fb, const XView& view, Range y) {
    fb.fill(theme_.bg);
    const Layout l = compute_layout(fb.width(), fb.height());
    const AxisTicks xt = make_ticks(view.x0, view.x1, xdivisions(l.plot.w));
    const AxisTicks yt = make_ticks(y.lo, y.hi, ydivisions(l.plot.h));

    // Axes (with grid) first so traces draw over the grid; the box is redrawn on top.
    draw_axes(fb, l, view.x0, view.x1, xt, y.lo, y.hi, yt, opts_.grid, theme_);
    std::vector<LegendEntry> legend;
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        const Signal& s = sigs_[i];
        const uint32_t color = s.color.value_or(theme_.trace_color(i));
        legend.push_back({s.name, color, s.style, s.visible});
        if (!s.visible) continue;
        const TraceStyle st{s.style, color, s.thickness > 0 ? s.thickness : opts_.thickness};
        paint_trace(fb, bins_[i], static_cast<int>(bins_[i].cols.size()), view, l.plot, y.lo, y.hi, st);
    }
    const Rect& p = l.plot;
    fb.rect_outline({p.x - 1, p.y - 1, p.w + 2, p.h + 2}, theme_.fg);
    if (opts_.legend && show_legend_) draw_legend(fb, p, legend, theme_);
    draw_title(fb, l, opts_.title, theme_);
    fb.text(Framebuffer::kCharW, l.readout.y + 4, cmode_name(cmode_), theme_.fg);
}

void TracePlot::render(Framebuffer& fb, const CancelToken& ct) {
    xview_ = autoscale_x();
    const XView view = make_view(xview_, compute_layout(fb.width(), fb.height()).plot.w);
    reduce(view, ct);
    yview_ = autoscale_y();
    paint(fb, view, yview_);
}

}  // namespace ssp
