/// @file content.h
/// @brief What a screen plots: the interface xplot (TraceContent) and xraster (RasterContent)
/// implement, plus the reduce request/result they share.
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "app/settings.h"
#include "core/cancel.h"
#include "core/component.h"
#include "render/frame.h"
#include "render/trace.h"

namespace ssp {

/// @brief Per-trace sample window [first, last) for IR mode.
using SampleRange = std::pair<std::size_t, std::size_t>;

/// @brief Everything a reduce depends on. Two equal requests give identical results.
struct ReduceRequest {
    XView view;                   ///< x0/x1 are ignored when auto_x is set.
    bool index = false;           ///< x = sample index instead of xstart + i*xdelta.
    bool auto_x = false;          ///< IR: autoscale x (the real part) from the samples.
    std::optional<Range> y;       ///< Fixed y range; empty = autoscale.
    std::vector<SampleRange> ir;  ///< IR: samples of each trace to plot (the time window).
    Reduce reduce = Reduce::Max;  ///< Raster: how a pixel combines its samples.
};

/// @brief Output of one reduce, written by the worker and read by the UI thread.
///
/// Everything except `done`, `ended`, `x`, `y` and the counters is sized before the job
/// starts. The worker writes results and then publishes a new `done` (release); the UI reads
/// only what is below the published count. `total` is what `done` counts up to: pixel
/// columns for traces, rows for rasters, a progress scale for IR.
struct ReduceResult {
    ReduceRequest req;
    XView view;                       ///< req.view with width/height clamped to >= 1.
    int total = 1;
    std::vector<std::size_t> traces;  ///< Trace indices reduced (those visible at submit time).
    std::vector<TraceBins> bins;      ///< Indexed by trace; empty for traces not reduced.
    std::vector<Range> xmap;          ///< Per trace {xstart, xdelta} used for this view.
    /// Ranges fixed for this result: from the zoom level, or computed by the worker before it
    /// publishes anything (dots need exact rows; IR autoscales both axes). Read once done > 0.
    std::optional<Range> x;
    std::optional<Range> y;
    /// IR only: per-trace hit counts, width*height cells, row-major from the top. Only the
    /// worker writes (relaxed load+store), so the UI can read counts while they grow.
    std::vector<std::unique_ptr<std::atomic<uint32_t>[]>> density;
    /// Raster only: reduced z per pixel, width*height, row-major from the top; NaN = no data.
    std::vector<float> zimg;
    std::atomic<int> done{0};
    std::atomic<bool> ended{false};   ///< Job exited, finished or cancelled.

    int width() const { return view.width; }
    bool complete() const { return done.load(std::memory_order_acquire) >= total; }
    double progress(int d) const { return static_cast<double>(d) / total; }
};

/// @brief What paint() needs besides the result.
struct PaintArgs {
    Rect plot;
    int done = 0;
    Range x, y, z;  ///< z: raster colour range.
    const Settings& set;
    const Theme& th;
};

class Content {
public:
    virtual ~Content() = default;

    virtual bool is_raster() const { return false; }
    /// @brief Settings to start with (from the options the content was made with).
    virtual Settings initial_settings() const = 0;
    /// @brief Level-0 view for `s`.
    virtual View home(const Settings& s) const = 0;
    virtual bool supports(CMode m) const { return m != CMode::Auto; }
    /// @brief {xstart, xdelta} of the first trace / raster columns, for index readouts.
    virtual Range xaxis() const = 0;

    /// @brief UI thread: allocate a result for `req`.
    virtual std::shared_ptr<ReduceResult> new_result(const ReduceRequest& req) const = 0;
    /// @brief True if `r` is complete and is exactly what new_result(req) would compute now.
    virtual bool reusable(const ReduceResult& r, const ReduceRequest& req) const = 0;
    /// @brief Worker: fill `r`, calling `progress` after publishing.
    /// @throws Cancelled when `ct` is cancelled.
    virtual void reduce(ReduceResult& r, const CancelToken& ct,
                        const std::function<void()>& progress = {}) const = 0;
    /// @brief Autoscale extent of what is finished: y for traces, z for rasters.
    virtual Span extent(const ReduceResult& r, int done) const = 0;
    virtual void paint(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const = 0;

    virtual std::vector<LegendEntry> legend(const Theme&) const { return {}; }
    /// @brief z (display units) of the sample at data point (x, y), rasters only.
    virtual std::optional<double> z_at(double, double, const Settings&) const { return std::nullopt; }
};

/// @brief Pad an autoscaled extent like SigPlot (2% each side); [-1, 1] when empty.
Range autoscale_y(const Span& extent);

}  // namespace ssp
