/// @file trace.h
/// @brief xplot content: reduce a trace to per-column bins (expensive), then paint them (cheap).
#pragma once

#include <atomic>
#include <cstddef>
#include <vector>

#include "core/cancel.h"
#include "core/component.h"
#include "core/lod.h"
#include "render/framebuffer.h"

namespace ssp {

/// @brief One pixel column: min/max plus the first and last sample (for connectors), in
/// plotted (display) units. count == 0 means no samples fall in the column.
struct Bin {
    double lo = 0.0, hi = 0.0, first = 0.0, last = 0.0;
    std::size_t count = 0;
};

/// @brief Point just outside the view, so lines run to the plot edge.
struct Context {
    bool valid = false;
    double x = 0.0, v = 0.0;
};

struct TraceBins {
    std::vector<Bin> cols;
    Context before, after;
    Span extent;  ///< y extent of all finished columns (display units).
};

/// @brief What a reduce depends on. A different XView means different bins.
struct XView {
    double x0 = 0.0, x1 = 1.0;  ///< Visible x range.
    int width = 1;              ///< Number of pixel columns.
    CMode cmode = CMode::Real;  ///< Resolved (never Auto or IR).
    PhaseUnits units = PhaseUnits::Radians;
};

/// @brief Fill `out` for samples of `lod` whose x = xstart + i*xdelta lies in the view.
/// Columns finish left to right; `done` (if given) is advanced after each one.
/// @throws Cancelled when `ct` is cancelled.
void reduce_trace(Lod& lod, double xstart, double xdelta, const XView& view, TraceBins& out,
                  const CancelToken& ct = {}, std::atomic<int>* done = nullptr);

struct TraceStyle {
    Style style = Style::Lines;
    uint32_t color = 0xFFFFFF;
    int thickness = 1;
};

/// @brief Draw the first `done` columns of `bins` into `plot`, mapping y in [y0, y1].
void paint_trace(Framebuffer& fb, const TraceBins& bins, int done, const XView& view, Rect plot,
                 double y0, double y1, const TraceStyle& style);

}  // namespace ssp
