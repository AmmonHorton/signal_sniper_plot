/// @file trace_content.h
/// @brief xplot data: signals, their pyramids, and reduce results that fill in over time.
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "app/content.h"
#include "app/settings.h"
#include "core/lod.h"
#include "render/frame.h"
#include "render/trace.h"
#include "ssp/plot.h"

namespace ssp {

class TraceContent final : public Content {
public:
    /// @brief IR mode draws sample-to-sample lines only up to this many samples in the window;
    /// beyond that only the density is drawn.
    static constexpr std::size_t kIrLineLimit = 50'000;

    /// @throws std::invalid_argument for unusable signals or options.
    TraceContent(std::vector<Signal> signals, const PlotOptions& options);

    Settings initial_settings() const override;
    View home(const Settings& s) const override;
    Range xaxis() const override { return {sigs_[0].xstart, sigs_[0].xdelta}; }
    /// @brief IR window to start with when the initial mode is IR: xrange or everything.
    std::vector<SampleRange> initial_ir() const;

    std::size_t size() const { return sigs_.size(); }
    const Signal& signal(std::size_t i) const { return sigs_[i]; }
    bool any_complex() const;

    /// @brief CMode to start in: options.cmode with Auto resolved.
    CMode initial_cmode() const { return initial_cmode_; }

    void set_visible(std::size_t i, bool v) { sigs_[i].visible = v; }
    /// @brief lines → dots → both → lines.
    void cycle_style(std::size_t i);

    /// @brief x range covering every visible trace (all traces if none are visible).
    Range x_extent(bool index) const;

    /// @brief Samples of each trace whose x lies in `x` (IR's time window).
    std::vector<SampleRange> samples_in(Range x, bool index) const;

    /// @brief Allocate for the traces visible now.
    std::shared_ptr<ReduceResult> new_result(const ReduceRequest& req) const override;
    /// @brief Also requires the same visible traces and dots styles.
    bool reusable(const ReduceResult& r, const ReduceRequest& req) const override;
    void reduce(ReduceResult& r, const CancelToken& ct,
                const std::function<void()>& progress = {}) const override;
    /// @brief y extent of columns [0, done) of the traces that are visible now.
    Span extent(const ReduceResult& r, int done) const override;
    void paint(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const override;
    /// @brief Legend rows (every trace), or none when the legend would say nothing useful.
    std::vector<LegendEntry> legend(const Theme& th) const override;

private:
    void reduce_ir(ReduceResult& r, const CancelToken& ct,
                   const std::function<void()>& progress) const;
    void paint_ir(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const;

    std::vector<Signal> sigs_;
    PlotOptions opts_;
    std::vector<std::unique_ptr<Lod>> lods_;  // used only by whichever thread runs reduce()
    CMode initial_cmode_ = CMode::Real;
    bool show_legend_ = false;
};

}  // namespace ssp
