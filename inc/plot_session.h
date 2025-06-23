// plot_session.h
#pragma once

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <vector>
#include <stack>
#include <optional>
#include <string>
#include "trace_utils.h"

namespace xplot {

enum class PlotMode { Magnitude, Real, Imag, Phase, IQ };
enum class XAxisMode { TIME, INDEX };

struct ZoomRegion {
    double xmin, xmax;
    double ymin, ymax;
};

struct Button {
    std::string label;
    int x, width;
};

struct LegendItem {
    int x, y, width, height;
};

struct DecimatedTrace {
    std::vector<std::pair<double, double>> y_minmax_per_pixel;
    ZoomRegion cached_region;
    PlotMode cached_mode;
    int cached_width = 0;
    bool dirty = true;

    void update_if_needed(const std::vector<PlotSample>& samples,
                                      PlotMode mode,
                                      const ZoomRegion& view,
                                      int pixel_width);
};

class PlotSession {
public:
    PlotSession(const std::string& title = "", int line_thickness = 1);

    void add_trace(const Trace& trace);
    void set_axis_mode(XAxisMode mode);
    void set_plot_mode(PlotMode mode);
    void set_zoom_range(ZoomRegion zr);
    void set_fixed_y_range(std::optional<std::pair<double, double>> ylimits);
    void set_fixed_x_range(std::optional<std::pair<double, double>> xlimits); // <-- add this
    
    void run();         // Interactive event loop
    void render_once(); // For headless or pre-rendered output
    
private:
    void reset_zoom();
    void init_x11();
    void handle_events();
    void render_pixmap();
    void draw_overlay();
    void draw_legend();
    void draw_toolbar();
    void draw_axes();
    void save_pixmap_to_ppm(const char* filename);

    ZoomRegion autoscale_region(const std::vector<PlotSample>&);

    // X11 state
    Display* dpy_ = nullptr;
    Window win_ = 0;
    GC gc_ = 0;
    Pixmap pixmap_ = 0;
    int screen_ = 0;
    int width_ = 1000;
    int height_ = 600;

    // Trace data and rendering state
    std::vector<Trace> traces_;
    std::vector<DecimatedTrace> decimated_cache_;
    std::vector<Button> toolbar_buttons_;
    std::vector<LegendItem> legend_boxes_;

    std::stack<ZoomRegion> zoom_stack_;
    ZoomRegion current_view_;
    int zoom_depth_ = 0;
    std::string readout_;

    PlotMode mode_ = PlotMode::Magnitude;
    XAxisMode xaxis_mode_ = XAxisMode::TIME;
    std::string title_;
    int line_thickness_ = 1;
    std::optional<std::pair<double, double>> fixed_y_range_;
    std::optional<std::pair<double, double>> fixed_x_range_;

    // Input / interaction
    bool dragging_ = false;
    int bx0_ = -1, by0_ = -1, bx1_ = -1, by1_ = -1;
    bool pixmap_dirty_ = true;
};

// Utility helpers (exposed for UI labels, etc.)
std::string mode_to_string(PlotMode mode);
double compute_value(const PlotSample& s, PlotMode mode);

}  // namespace xplot
