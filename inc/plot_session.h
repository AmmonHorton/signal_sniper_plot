/// @file plot_session.h
/// @brief PlotSession and supporting types for the interactive X11 plot window.
#pragma once

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <vector>
#include <stack>
#include <optional>
#include <string>
#include "trace_utils.h"

namespace xplot {

/// @brief Determines which signal component is rendered on the Y-axis.
enum class PlotMode {
    Magnitude, ///< sqrt(real² + imag²).
    Real,      ///< Real component only.
    Imag,      ///< Imaginary component only.
    Phase,     ///< atan2(real, imag) in radians.
    IQ         ///< Scatter plot of imag vs. real.
};

/// @brief Controls how x-axis tick labels are formatted.
enum class XAxisMode {
    TIME,  ///< Labels show the raw x value of each tick.
    INDEX  ///< Labels show the nearest integer sample index.
};

/// @brief Axis-aligned rectangle defining the visible data range.
struct ZoomRegion {
    double xmin; ///< Left x boundary.
    double xmax; ///< Right x boundary.
    double ymin; ///< Bottom y boundary.
    double ymax; ///< Top y boundary.
};

/// @brief Toolbar button geometry and label, populated by draw_toolbar().
struct Button {
    std::string label; ///< Text shown on the button face.
    int x;             ///< Left pixel edge of the button.
    int width;         ///< Pixel width of the button.
};

/// @brief Bounding box of a legend colour swatch, used for hit-testing clicks.
struct LegendItem {
    int x;      ///< Left pixel edge.
    int y;      ///< Top pixel edge.
    int width;  ///< Width of the clickable area.
    int height; ///< Height of the clickable area.
};

/// @brief Pixel geometry of the plot's sub-regions, recomputed on each resize.
struct PlotLayout {
    int px = 0, py = 0;    ///< Data area top-left corner (pixels).
    int pw = 1, ph = 1;    ///< Data area width/height; clamped to >= 1.
    int legend_y  = 0;     ///< Top edge of the legend row (pixels).
    int toolbar_y = 0;     ///< Top edge of the toolbar row (pixels).
    int readout_x = 0;     ///< Left edge of the mode/cursor readout text.
};

/// @brief Per-trace state carried across streaming render frames.
/// Tracks the last drawn pixel position and the next sample index to consume.
struct TraceRenderState {
    int    prev_x        = -1;   ///< X pixel of the last drawn point; -1 if none yet.
    int    prev_y        = -1;   ///< Y pixel of the last drawn point; -1 if none yet.
    double prev_data_y   = 0.0;  ///< Data-space Y value at prev_x / prev_y.
    int    last_oor_x    = -1;   ///< X pixel of the last out-of-range bin; -1 if none.
    double last_oor_val  = 0.0;  ///< Extremal data value of the last out-of-range bin.
    size_t sample_cursor = 0;    ///< Index of the next sample to consume during streaming.
};

/// @brief Per-trace decimation cache mapping each pixel column to a [min, max] value range.
/// Marked dirty when the view, mode, or pixel width changes; rebuilt incrementally during streaming.
struct DecimatedTrace {
    std::vector<std::pair<double, double>> y_minmax_per_pixel; ///< Per-column [min, max] pairs; length == pixel_width.
    ZoomRegion cached_region; ///< View region the cache was built for.
    PlotMode   cached_mode;   ///< Plot mode the cache was built for.
    int  cached_width = 0;    ///< Pixel width the cache was built for.
    bool dirty        = true; ///< True if the cache must be rebuilt before use.

    /// @brief Populate the cache from @p samples for the given view and pixel width.
    /// A no-op when the cache is already valid (not dirty and all parameters unchanged).
    void update_if_needed(const std::vector<PlotSample>& samples,
                          PlotMode mode,
                          const ZoomRegion& view,
                          int pixel_width);
};

/// @brief Interactive X11 plotting session.
/// Manages window lifecycle, event handling, zoom/pan, and incremental per-trace streaming.
class PlotSession {
public:
    /// @brief Construct a session with an optional title and line thickness.
    PlotSession(const std::string& title = "", int line_thickness = 1);

    /// @brief Add a trace to be plotted. The trace is copied into the session.
    void add_trace(const Trace& trace);

    /// @brief Set the x-axis label mode (time values or sample indices).
    void set_axis_mode(XAxisMode mode);

    /// @brief Set the initial plot mode (which signal component to display).
    void set_plot_mode(PlotMode mode);

    /// @brief Override the initial view region instead of auto-scaling from data.
    void set_zoom_range(ZoomRegion zr);

    /// @brief Fix the Y-axis range; pass nullopt to re-enable auto-scaling.
    void set_fixed_y_range(std::optional<std::pair<double, double>> ylimits);

    /// @brief Fix the X-axis range; pass nullopt to re-enable auto-scaling.
    void set_fixed_x_range(std::optional<std::pair<double, double>> xlimits);

    /// @brief Open the window and enter the interactive event loop. Blocks until the window is closed.
    void run();

    /// @brief Render all traces into the pixmap once and flush without entering the event loop.
    void render_once();

private:
    /// @brief Reset the zoom stack and re-autoscale from the first trace.
    void reset_zoom();

    /// @brief Open the X11 display, create the window, GC, and backing pixmap.
    void init_x11();

    /// @brief Clear the pixmap, invalidate stale caches, draw static UI elements, and seed trace streaming.
    void render_pixmap_init();

    /// @brief Render up to @p n pixel columns of the currently active streaming trace.
    /// Uses the decimation cache when valid; otherwise populates it column-by-column from samples.
    void render_stream_columns(int n);

    /// @brief Initialize streaming state for trace @p t.
    /// Resets TraceRenderState, seeds sample_cursor via binary search, and allocates bins when dirty.
    void start_trace_stream(size_t t);

    /// @brief Draw legend colour swatches and labels into the pixmap.
    void draw_legend(const PlotLayout& l);

    /// @brief Draw toolbar button outlines and labels into the pixmap.
    void draw_toolbar(const PlotLayout& l);

    /// @brief Draw x- and y-axis tick marks and numeric labels into the pixmap.
    void draw_axes(const PlotLayout& l);

    /// @brief Save the current backing pixmap as a raw PPM file.
    void save_pixmap_to_ppm(const char* filename);

    /// @brief Compute pixel geometry for the data area, legend, and toolbar from the current window size.
    PlotLayout compute_layout() const;

    /// @brief Compute an auto-scaled ZoomRegion from the given sample set.
    /// Pads each axis by 5%; respects fixed_x_range_ / fixed_y_range_ when set.
    ZoomRegion autoscale_region(const std::vector<PlotSample>&);

    // ── X11 state ─────────────────────────────────────────────────────────────
    Display* dpy_    = nullptr; ///< Connection to the X11 display server.
    Window   win_    = 0;       ///< Top-level application window.
    GC       gc_     = 0;       ///< Graphics context used for all drawing.
    Pixmap   pixmap_ = 0;       ///< Off-screen backing store; blitted to win_ each frame.
    Atom     wm_delete_window_ = 0; ///< WM_DELETE_WINDOW atom for clean close handling.
    int      screen_ = 0;       ///< Default screen index.
    int      width_  = 1000;    ///< Current window width in pixels.
    int      height_ = 600;     ///< Current window height in pixels.

    // ── Trace data and rendering state ────────────────────────────────────────
    std::vector<Trace>          traces_;          ///< All traces added to this session.
    std::vector<DecimatedTrace> decimated_cache_; ///< One cache entry per trace.
    std::vector<Button>         toolbar_buttons_; ///< Populated by draw_toolbar().
    std::vector<LegendItem>     legend_boxes_;    ///< Populated by draw_legend(); used for click hit-testing.

    PlotLayout layout_; ///< Current pixel layout; recomputed on resize and each render.

    std::stack<ZoomRegion> zoom_stack_;  ///< History of zoom regions for right-click undo.
    ZoomRegion current_view_;            ///< Active view region in data coordinates.
    int        zoom_depth_ = 0;          ///< Current zoom stack depth; capped at MAX_ZOOM_HISTORY.
    std::string readout_;                ///< Cursor position string rendered in the toolbar area.

    PlotMode  mode_       = PlotMode::Magnitude; ///< Active plot mode.
    XAxisMode xaxis_mode_ = XAxisMode::TIME;     ///< Active x-axis label mode.
    std::string title_;       ///< Window title drawn above the data area.
    int line_thickness_ = 1;  ///< Pixel thickness for LINES-style traces.

    std::optional<std::pair<double, double>> fixed_y_range_; ///< User-supplied Y bounds; empty → auto-scale.
    std::optional<std::pair<double, double>> fixed_x_range_; ///< User-supplied X bounds; empty → auto-scale.

    // ── Streaming render state ─────────────────────────────────────────────────
    int    stream_xpix_         = 0;     ///< Next pixel column to render for the active trace.
    bool   streaming_active_    = false; ///< True while at least one trace is still being streamed.
    bool   stream_paused_       = false; ///< True when streaming is paused via middle-click.
    size_t stream_active_trace_ = 0;    ///< Index of the trace currently being streamed.
    std::vector<size_t>           stream_pending_;     ///< FIFO of trace indices queued after the active one.
    std::vector<TraceRenderState> stream_trace_state_; ///< Per-trace streaming state; size == traces_.size().

    // ── Input / interaction ────────────────────────────────────────────────────
    bool dragging_ = false;             ///< True while the user holds the left button for a zoom drag.
    int bx0_ = -1, by0_ = -1;          ///< Zoom-drag start pixel.
    int bx1_ = -1, by1_ = -1;          ///< Zoom-drag current/end pixel.
    int mouse_x_ = -1, mouse_y_ = -1;  ///< Most recent cursor position for crosshair rendering.
    bool pixmap_dirty_ = true;          ///< True when the pixmap must be fully re-rendered before the next blit.
};

/// @brief Convert a PlotMode enumerator to its short display string (e.g., PlotMode::IQ → "I/Q").
std::string mode_to_string(PlotMode mode);

/// @brief Extract the scalar Y value from a PlotSample for the given PlotMode.
/// Returns 0 for IQ mode (scatter plots are handled separately).
double compute_value(const PlotSample& s, PlotMode mode);

}  // namespace xplot
