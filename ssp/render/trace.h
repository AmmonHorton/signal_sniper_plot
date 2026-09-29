/// @file trace.h
/// @brief xplot content: reduce a trace to per-column bins (expensive), then paint them (cheap).
#pragma once

#include <atomic>
#include <cmath>
#include <cstddef>
#include <vector>

#include "ssp/core/cancel.h"
#include "ssp/core/component.h"
#include "ssp/core/lod.h"
#include "ssp/render/framebuffer.h"

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

    /// Dots mode only: which pixel rows of each column hold at least one sample, so dots are
    /// drawn exactly where samples are (never filled in between). Bit r of column c is
    /// occ[c * occ_words + r / 64] >> (r % 64). Empty when not requested; rows are for y
    /// range `occ_y`, which the reducer's owner sets before any column is published.
    std::vector<uint64_t> occ;
    int occ_words = 0;
    Range occ_y;

    bool has_row(int c, int r) const { return occ[c * occ_words + r / 64] >> (r % 64) & 1; }
};

/// @brief Pixel row (0 = top) of value v in a data area `height` rows tall showing [y0, y1];
/// -1 above, `height` below.
inline int row_of(double v, double y0, double y1, int height) {
    const double r = std::floor((y1 - v) / (y1 - y0) * height);
    if (!(r >= 0.0)) return r < 0.0 ? -1 : height;  // NaN counts as off-screen
    return r >= height ? height : static_cast<int>(r);
}

/// @brief y extent (display units) of columns [0, done).
Span bins_extent(const TraceBins& bins, int done);

/// @brief What a reduce depends on. A different XView means different bins.
struct XView {
    double x0 = 0.0, x1 = 1.0;  ///< Visible x range.
    int width = 1;              ///< Number of pixel columns.
    int height = 1;             ///< Data-area rows (only dots-mode occupancy depends on it).
    CMode cmode = CMode::Real;  ///< Resolved (never Auto or IR).
    PhaseUnits units = PhaseUnits::Radians;
};

/// @brief Incrementally fills `out` for samples of `lod` whose x = xstart + i*xdelta lies in
/// the view, one column at a time from the left.
///
/// `out.cols` must already hold view.width default Bins. Only elements are written (never the
/// vector itself), so another thread may read columns below a count published after
/// run_to() returns. `before` is set by the constructor and `after` by the run_to() call that
/// reaches the last column.
class TraceReducer {
public:
    TraceReducer(Lod& lod, double xstart, double xdelta, const XView& view, TraceBins& out);

    /// @throws Cancelled when `ct` is cancelled.
    void run_to(int col_end, const CancelToken& ct = {});

    /// @brief y extent (display units) of every sample in the view. Cheap once the pyramid
    /// exists; used to fix the y range up front when dots need exact rows.
    Span view_extent(const CancelToken& ct = {}) const;
    int next() const { return next_; }

private:
    std::size_t boundary(int c) const;
    double value(std::size_t i) const;
    void mark_rows(int c, std::size_t i0, std::size_t i1, const Bin& b, const CancelToken& ct);

    Lod& lod_;
    double xstart_, xdelta_;
    XView view_;
    TraceBins& out_;
    Comp comp_;
    int next_ = 0;
    std::size_t prev_ = 0;
};

/// @brief Synchronous convenience: size `out` and reduce every column.
void reduce_trace(Lod& lod, double xstart, double xdelta, const XView& view, TraceBins& out,
                  const CancelToken& ct = {});

struct TraceStyle {
    Style style = Style::Lines;
    uint32_t color = 0xFFFFFF;
    int thickness = 1;
};

/// @brief Draw the first `done` columns of `bins` into `plot`, mapping y in [y0, y1].
void paint_trace(Framebuffer& fb, const TraceBins& bins, int done, const XView& view, Rect plot,
                 double y0, double y1, const TraceStyle& style);

}  // namespace ssp
