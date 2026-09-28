/// @file trace_content.h
/// @brief xplot data: signals, their pyramids, and reduce results that fill in over time.
#pragma once

#include <atomic>
#include <functional>
#include <memory>
#include <vector>

#include "app/settings.h"
#include "core/lod.h"
#include "render/frame.h"
#include "render/trace.h"
#include "ssp/plot.h"

namespace ssp {

/// @brief Bins for one view, written by the worker and read by the UI thread.
///
/// Everything except `done` and `ended` is sized before the job starts. The worker writes
/// columns [done, new_done) of each trace and then publishes new_done; the UI reads only
/// columns below the published count.
struct ReduceResult {
    XView view;
    std::vector<std::size_t> traces;  ///< Trace indices reduced (those visible at submit time).
    std::vector<TraceBins> bins;      ///< Indexed by trace; empty for traces not reduced.
    std::vector<Range> xmap;          ///< Per trace {xstart, xdelta} used for this view.
    std::optional<Range> y_requested;  ///< The zoom level's fixed y (reuse key), as submitted.
    /// y range fixed for this result: from the zoom level, or (when a dots trace needs exact
    /// pixel rows) computed by the worker before it publishes the first column. Read it only
    /// once done > 0.
    std::optional<Range> y;
    std::atomic<int> done{0};         ///< Finished columns (all traces).
    std::atomic<bool> ended{false};   ///< Job exited, finished or cancelled.

    int width() const { return view.width; }
    bool complete() const { return done.load(std::memory_order_acquire) >= view.width; }
};

class TraceContent {
public:
    /// @throws std::invalid_argument for unusable signals or options.
    TraceContent(std::vector<Signal> signals, const PlotOptions& options);

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

    /// @brief UI thread: allocate a result for `view` covering the traces visible now.
    /// @param y the zoom level's fixed y range, if it has one.
    std::shared_ptr<ReduceResult> new_result(const XView& view, bool index,
                                             std::optional<Range> y = std::nullopt) const;

    /// @brief True if `r` is complete and is exactly what new_result(view, index, y) would
    /// compute now (same view, mode, axes, y, visible traces and dots styles).
    bool reusable(const ReduceResult& r, const XView& view, bool index,
                  std::optional<Range> y) const;

    /// @brief Worker: fill `r` left to right, calling `progress` after publishing columns.
    /// @throws Cancelled when `ct` is cancelled.
    void reduce(ReduceResult& r, const CancelToken& ct,
                const std::function<void()>& progress = {}) const;

    /// @brief y extent of columns [0, done) of the traces that are visible now.
    Span extent(const ReduceResult& r, int done) const;

    void paint(Framebuffer& fb, const Rect& plot, const ReduceResult& r, int done, Range y,
               int thickness, const Theme& th) const;

    /// @brief Legend rows (every trace), or none when the legend would say nothing useful.
    std::vector<LegendEntry> legend(const Theme& th) const;

private:
    std::vector<Signal> sigs_;
    std::vector<std::unique_ptr<Lod>> lods_;  // used only by whichever thread runs reduce()
    CMode initial_cmode_ = CMode::Real;
    bool show_legend_ = false;
};

/// @brief Pad an autoscaled y extent like SigPlot (2% each side); [-1, 1] when empty.
Range autoscale_y(const Span& extent);

}  // namespace ssp
