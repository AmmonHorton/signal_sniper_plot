// IR (imaginary vs real) mode: a hit-count density of the samples in the chosen time window.

#include <chrono>
#include <cmath>

#include "app/trace_content.h"
#include "core/scan.h"

namespace ssp {
namespace {

constexpr std::size_t kIrChunk = 1 << 18;  // samples between progress checks
constexpr auto kPublishEvery = std::chrono::milliseconds(8);
constexpr double kMinBrightness = 0.35;    // a single hit is still clearly visible

uint32_t scale_color(uint32_t c, double f) {
    const auto ch = [&](int shift) {
        return static_cast<uint32_t>(std::lround(((c >> shift) & 0xff) * f)) << shift;
    };
    return ch(16) | ch(8) | ch(0);
}

}  // namespace

void TraceContent::reduce_ir(ReduceResult& r, const CancelToken& ct,
                             const std::function<void()>& progress) const {
    auto window = [&](std::size_t t) {
        return t < r.req.ir.size() ? r.req.ir[t] : SampleRange{0, sigs_[t].n};
    };
    if (!r.x || !r.y) {  // autoscale from the window's samples: two pyramid queries per trace
        Span re, im;
        for (std::size_t t : r.traces) {
            const auto [i0, i1] = window(t);
            re.merge(lods_[t]->query(Comp::Re, i0, i1, ct));
            im.merge(lods_[t]->query(Comp::Im, i0, i1, ct));
        }
        if (!r.x) r.x = autoscale_y(re);
        if (!r.y) r.y = autoscale_y(im);
    }

    const int w = r.view.width, h = r.view.height;
    const Range x = *r.x, y = *r.y;
    const double sx = w / (x.hi - x.lo), sy = h / (y.hi - y.lo);
    std::size_t total = 0, processed = 0;
    for (std::size_t t : r.traces) total += window(t).second - window(t).first;
    auto last_publish = std::chrono::steady_clock::now();

    for (std::size_t t : r.traces) {
        std::atomic<uint32_t>* d = r.density[t].get();
        const auto [i0, i1] = window(t);
        for (std::size_t c = i0; c < i1; c += kIrChunk) {
            const std::size_t e = std::min(i1, c + kIrChunk);
            for_each_sample(sigs_[t], c, e, ct, [&](double re, double im) {
                const double col = std::floor((re - x.lo) * sx), row = std::floor((y.hi - im) * sy);
                if (col >= 0.0 && col < w && row >= 0.0 && row < h) {
                    auto& cell = d[static_cast<std::size_t>(row) * w + static_cast<std::size_t>(col)];
                    cell.store(cell.load(std::memory_order_relaxed) + 1, std::memory_order_relaxed);
                }
            });
            processed += e - c;
            const auto now = std::chrono::steady_clock::now();
            if (now - last_publish >= kPublishEvery) {
                // Progress as a fraction of the width; never "complete" until the end.
                const int d_cols = static_cast<int>(double(w) * processed / std::max<std::size_t>(total, 1));
                r.done.store(std::clamp(d_cols, 1, w - 1), std::memory_order_release);
                last_publish = now;
                if (progress) progress();
            }
        }
    }
    r.done.store(w, std::memory_order_release);
    if (progress) progress();
}

void TraceContent::paint_ir(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const {
    const Rect& plot = a.plot;
    const Range x = a.x, y = a.y;
    const Theme& th = a.th;
    const int thickness = a.set.thickness;
    if (a.done <= 0 || !(x.hi > x.lo) || !(y.hi > y.lo)) return;
    const int w = r.view.width, h = r.view.height;
    const bool grid_matches = w == plot.w && h == plot.h;
    for (std::size_t t : r.traces) {
        const Signal& s = sigs_[t];
        if (!s.visible) continue;
        const uint32_t color = s.color.value_or(th.trace_color(t));
        const SampleRange win = t < r.req.ir.size() ? r.req.ir[t] : SampleRange{0, s.n};
        const bool few = win.second - win.first <= kIrLineLimit;
        const bool lines = s.style != Style::Dots && few;
        const bool density = s.style != Style::Lines || !few;

        if (density && grid_matches) {
            const std::atomic<uint32_t>* d = r.density[t].get();
            const std::size_t cells = static_cast<std::size_t>(w) * h;
            uint32_t most = 0;
            for (std::size_t k = 0; k < cells; ++k) most = std::max(most, d[k].load(std::memory_order_relaxed));
            const double norm = most > 1 ? 1.0 / std::log1p(double(most)) : 0.0;
            for (int row = 0; row < h; ++row) {
                uint32_t* px = fb.row(plot.y + row) + plot.x;
                for (int col = 0; col < w; ++col) {
                    const uint32_t n = d[static_cast<std::size_t>(row) * w + col].load(std::memory_order_relaxed);
                    if (n == 0) continue;
                    const double f = kMinBrightness + (1.0 - kMinBrightness) * std::log1p(double(n)) * norm;
                    px[col] = scale_color(color, most > 1 ? f : 1.0);
                }
            }
        }
        if (lines) {  // the actual sample-to-sample path, for small windows
            const int thick = s.thickness > 0 ? s.thickness : thickness;
            auto to_px = [&](std::size_t i, double& px, double& py) {
                px = plot.x + (sample_value(s, Comp::Re, i) - x.lo) / (x.hi - x.lo) * plot.w;
                py = plot.y + (y.hi - sample_value(s, Comp::Im, i)) / (y.hi - y.lo) * plot.h;
            };
            double px0 = 0, py0 = 0, px1, py1;
            for (std::size_t i = win.first; i < win.second; ++i) {
                to_px(i, px1, py1);
                if (i > win.first) fb.line(px0, py0, px1, py1, color, thick, plot);
                px0 = px1;
                py0 = py1;
            }
        }
    }
}

}  // namespace ssp
