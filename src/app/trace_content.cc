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
    : sigs_(std::move(signals)), opts_(o) {
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
}

Settings TraceContent::initial_settings() const {
    Settings s;
    s.title = opts_.title;
    s.cmode = initial_cmode_;
    s.phunits = opts_.phunits;
    s.index = opts_.index;
    s.grid = opts_.grid;
    s.legend = opts_.legend;
    s.thickness = opts_.thickness;
    return s;
}

View TraceContent::home(const Settings& s) const {
    View v;
    if (s.cmode == CMode::IR) {  // both axes autoscale from the window's samples
        v.auto_x = true;
        return v;
    }
    const bool original_axes = s.index == opts_.index;
    v.x = (original_axes && opts_.xrange) ? *opts_.xrange : x_extent(s.index);
    if (opts_.yrange && s.cmode == initial_cmode_) v.y = *opts_.yrange;
    return v;
}

std::vector<SampleRange> TraceContent::initial_ir() const {
    return samples_in(opts_.xrange.value_or(x_extent(opts_.index)), opts_.index);
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

std::vector<SampleRange> TraceContent::samples_in(Range x, bool index) const {
    std::vector<SampleRange> out;
    for (const Signal& sg : sigs_) {
        const double x0 = index ? 0.0 : sg.xstart, dx = index ? 1.0 : sg.xdelta;
        auto clamp_index = [&](double t) -> std::size_t {
            if (!(t > 0.0)) return 0;
            return t >= static_cast<double>(sg.n) ? sg.n : static_cast<std::size_t>(t);
        };
        const std::size_t i0 = clamp_index(std::ceil((x.lo - x0) / dx));
        const std::size_t i1 = clamp_index(std::floor((x.hi - x0) / dx) + 1.0);
        out.push_back({i0, std::max(i0, i1)});
    }
    return out;
}

std::shared_ptr<ReduceResult> TraceContent::new_result(const ReduceRequest& req) const {
    auto r = std::make_shared<ReduceResult>();
    r->req = req;
    r->view = req.view;
    r->view.width = std::max(1, req.view.width);
    r->view.height = std::max(1, req.view.height);
    r->total = r->view.width;
    if (!req.auto_x) r->x = Range{req.view.x0, req.view.x1};
    r->y = req.y;
    r->bins.resize(sigs_.size());
    r->xmap.resize(sigs_.size());
    const bool ir = req.view.cmode == CMode::IR;
    if (ir) r->density.resize(sigs_.size());
    const std::size_t cells = static_cast<std::size_t>(r->view.width) * r->view.height;
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        r->xmap[i] = req.index ? Range{0.0, 1.0} : Range{sigs_[i].xstart, sigs_[i].xdelta};
        if (!sigs_[i].visible) continue;
        r->traces.push_back(i);
        if (ir) {
            r->density[i].reset(new std::atomic<uint32_t>[cells]);
            for (std::size_t k = 0; k < cells; ++k) r->density[i][k].store(0, std::memory_order_relaxed);
            continue;
        }
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

bool TraceContent::reusable(const ReduceResult& r, const ReduceRequest& req) const {
    if (!r.complete()) return false;
    const ReduceRequest& q = r.req;
    const XView &a = q.view, &b = req.view;
    const bool ir = b.cmode == CMode::IR;
    if (a.width != b.width || a.height != b.height || a.cmode != b.cmode || a.units != b.units ||
        q.index != req.index || q.auto_x != req.auto_x || !same(q.y, req.y)) {
        return false;
    }
    if (!req.auto_x && (a.x0 != b.x0 || a.x1 != b.x1)) return false;
    if (ir && q.ir != req.ir) return false;
    std::size_t k = 0;  // r.traces must be exactly the visible traces, with matching styles
    for (std::size_t i = 0; i < sigs_.size(); ++i) {
        if (!sigs_[i].visible) continue;
        if (k >= r.traces.size() || r.traces[k++] != i) return false;
        if (!ir && r.bins[i].occ.empty() != (sigs_[i].style != Style::Dots)) return false;
    }
    return k == r.traces.size();
}

void TraceContent::reduce(ReduceResult& r, const CancelToken& ct,
                          const std::function<void()>& progress) const {
    if (r.view.cmode == CMode::IR) return reduce_ir(r, ct, progress);
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

void TraceContent::paint(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const {
    if (r.view.cmode == CMode::IR) return paint_ir(fb, r, a);
    for (std::size_t t : r.traces) {
        const Signal& s = sigs_[t];
        if (!s.visible) continue;
        const TraceStyle st{s.style, s.color.value_or(a.th.trace_color(t)),
                            s.thickness > 0 ? s.thickness : a.set.thickness};
        paint_trace(fb, r.bins[t], a.done, r.view, a.plot, a.y.lo, a.y.hi, st);
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
