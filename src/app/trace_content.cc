#include "app/trace_content.h"

#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ssp {
namespace {

constexpr double kYPad = 0.02;           // SigPlot pads autoscaled y by 2% each side
constexpr int kColumnsPerStep = 8;       // columns per trace between progress checks
constexpr auto kPublishEvery = std::chrono::milliseconds(8);

void check_range(const std::optional<Range>& r, const char* what) {
    if (r && !(std::isfinite(r->lo) && std::isfinite(r->hi) && r->hi > r->lo)) {
        throw std::invalid_argument(std::string(what) + " must be finite with hi > lo");
    }
}

}  // namespace

Range autoscale_y(const Span& extent) {
    if (extent.empty()) return {-1.0, 1.0};
    if (!(extent.hi > extent.lo)) return {extent.lo - 1.0, extent.hi + 1.0};
    const double pad = kYPad * (extent.hi - extent.lo);
    return {extent.lo - pad, extent.hi + pad};
}

TraceContent::TraceContent(std::vector<Signal> signals, const PlotOptions& o)
    : sigs_(std::move(signals)) {
    if (sigs_.empty()) throw std::invalid_argument("plot needs at least one signal");
    check_range(o.xrange, "xrange");
    check_range(o.yrange, "yrange");
    if (o.thickness < 1) throw std::invalid_argument("thickness must be >= 1");

    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        Signal& s = sigs_[i];
        const std::string who = "signal " + std::to_string(i) + ": ";
        if (s.n > 0 && !s.data) throw std::invalid_argument(who + "data is null");
        if (s.stride < 1) throw std::invalid_argument(who + "stride must be >= 1");
        if (!(std::isfinite(s.xdelta) && s.xdelta > 0.0)) {
            throw std::invalid_argument(who + "xdelta must be finite and > 0");
        }
        if (!std::isfinite(s.xstart)) throw std::invalid_argument(who + "xstart must be finite");
        show_legend_ = show_legend_ || !s.name.empty() || sigs_.size() > 1;
        if (s.name.empty()) s.name = "Trace " + std::to_string(i);
        lods_.push_back(std::make_unique<Lod>(s));
    }

    initial_cmode_ = o.cmode;
    if (initial_cmode_ == CMode::Auto) initial_cmode_ = any_complex() ? CMode::Mag : CMode::Real;
    if (initial_cmode_ == CMode::IR) {
        throw std::invalid_argument("IR (imag vs real) mode is not implemented yet");
    }
}

bool TraceContent::any_complex() const {
    for (const auto& s : sigs_) {
        if (s.complex) return true;
    }
    return false;
}

void TraceContent::cycle_style(std::size_t i) {
    Style& st = sigs_[i].style;
    st = st == Style::Lines ? Style::Dots : st == Style::Dots ? Style::LinesDots : Style::Lines;
}

Range TraceContent::x_extent(bool index) const {
    Span s;
    bool any_visible = false;
    for (const auto& sg : sigs_) any_visible = any_visible || (sg.visible && sg.n > 0);
    for (const auto& sg : sigs_) {
        if (sg.n == 0 || (any_visible && !sg.visible)) continue;
        const double x0 = index ? 0.0 : sg.xstart, dx = index ? 1.0 : sg.xdelta;
        s.add(x0);
        s.add(x0 + (sg.n - 1) * dx);
    }
    if (s.empty()) return {-1.0, 1.0};
    if (!(s.hi > s.lo)) return {s.lo - 1.0, s.hi + 1.0};
    return {s.lo, s.hi};
}

std::shared_ptr<ReduceResult> TraceContent::new_result(const XView& view, bool index,
                                                     std::optional<Range> y) const {
    auto r = std::make_shared<ReduceResult>();
    r->view = view;
    r->view.width = std::max(1, view.width);
    r->view.height = std::max(1, view.height);
    r->y = y;
    r->y_requested = y;
    r->bins.resize(sigs_.size());
    r->xmap.resize(sigs_.size());
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        r->xmap[i] = index ? Range{0.0, 1.0} : Range{sigs_[i].xstart, sigs_[i].xdelta};
        if (!sigs_[i].visible) continue;
        r->traces.push_back(i);
        TraceBins& b = r->bins[i];
        b.cols.assign(r->view.width, Bin{});
        if (sigs_[i].style == Style::Dots) {
            b.occ_words = (r->view.height + 63) / 64;
            b.occ.assign(static_cast<std::size_t>(r->view.width) * b.occ_words, 0);
        }
    }
    return r;
}

namespace {
bool same(const std::optional<Range>& a, const std::optional<Range>& b) {
    if (a.has_value() != b.has_value()) return false;
    return !a || (a->lo == b->lo && a->hi == b->hi);
}
}  // namespace

bool TraceContent::reusable(const ReduceResult& r, const XView& view, bool index,
                            std::optional<Range> y) const {
    if (!r.complete()) return false;
    const XView& v = r.view;
    if (v.x0 != view.x0 || v.x1 != view.x1 || v.width != std::max(1, view.width) ||
        v.height != std::max(1, view.height) || v.cmode != view.cmode || v.units != view.units ||
        !same(r.y_requested, y)) {
        return false;
    }
    std::size_t k = 0;  // r.traces must be exactly the visible traces, with matching styles
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        if (!sigs_[i].visible) continue;
        if (k >= r.traces.size() || r.traces[k++] != i) return false;
        if (r.bins[i].occ.empty() != (sigs_[i].style != Style::Dots)) return false;
        const Range xm = index ? Range{0.0, 1.0} : Range{sigs_[i].xstart, sigs_[i].xdelta};
        if (r.xmap[i].lo != xm.lo || r.xmap[i].hi != xm.hi) return false;
    }
    return k == r.traces.size();
}

void TraceContent::reduce(ReduceResult& r, const CancelToken& ct,
                          const std::function<void()>& progress) const {
    std::vector<TraceReducer> reducers;
    for (std::size_t t : r.traces) {
        reducers.emplace_back(*lods_[t], r.xmap[t].lo, r.xmap[t].hi, r.view, r.bins[t]);
    }
    // Exact dot rows need the final y range before any column is marked. With a fixed range
    // that is free; otherwise autoscale from the whole view first (one pyramid query per trace).
    bool dots = false;
    for (std::size_t t : r.traces) dots = dots || !r.bins[t].occ.empty();
    if (dots) {
        if (!r.y) {
            Span all;
            for (const auto& red : reducers) all.merge(red.view_extent(ct));
            r.y = autoscale_y(all);
        }
        for (std::size_t t : r.traces) r.bins[t].occ_y = *r.y;
    }
    // All traces advance together so the whole plot sweeps left to right.
    auto last_publish = std::chrono::steady_clock::now();
    const int w = r.width();
    for (int c = 0; c < w;) {
        const int end = std::min(w, c + kColumnsPerStep);
        for (auto& red : reducers) red.run_to(end, ct);
        c = end;
        const auto now = std::chrono::steady_clock::now();
        if (c == w || now - last_publish >= kPublishEvery) {
            r.done.store(c, std::memory_order_release);
            last_publish = now;
            if (progress) progress();
        }
    }
}

Span TraceContent::extent(const ReduceResult& r, int done) const {
    Span s;
    for (std::size_t t : r.traces) {
        if (sigs_[t].visible) s.merge(bins_extent(r.bins[t], done));
    }
    return s;
}

void TraceContent::paint(Framebuffer& fb, const Rect& plot, const ReduceResult& r, int done,
                         Range y, int thickness, const Theme& th) const {
    for (std::size_t t : r.traces) {
        const Signal& s = sigs_[t];
        if (!s.visible) continue;
        const TraceStyle st{s.style, s.color.value_or(th.trace_color(t)),
                            s.thickness > 0 ? s.thickness : thickness};
        paint_trace(fb, r.bins[t], done, r.view, plot, y.lo, y.hi, st);
    }
}

std::vector<LegendEntry> TraceContent::legend(const Theme& th) const {
    std::vector<LegendEntry> out;
    if (!show_legend_) return out;
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        const Signal& s = sigs_[i];
        out.push_back({s.name, s.color.value_or(th.trace_color(i)), s.style, s.visible});
    }
    return out;
}

}  // namespace ssp
