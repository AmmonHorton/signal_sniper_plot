#include "plot.hpp"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <stack>
#include <sstream>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <algorithm>
#include <climits>
#include <map>

namespace xplot {

enum class PlotMode { Magnitude, Real, Imag, Phase };
enum class XAxisMode { TIME, INDEX };
enum class PlotRenderStyle { LINES, DOTS };

struct PlotSample {
    double time;
    double real;
    double imag;
};

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

static std::vector<PlotSample> g_samples;
static std::stack<ZoomRegion> g_zoom_stack;
static int g_zoom_depth = 0;
constexpr int MAX_ZOOM_HISTORY = 5;

static std::vector<PlotRenderStyle> g_trace_styles;
static std::vector<bool> g_trace_visibility;
static std::vector<LegendItem> g_legend_boxes;

static std::vector<std::vector<PlotSample>> g_traces;
static PlotMode g_mode = PlotMode::Magnitude;
static XAxisMode g_xaxis_mode = XAxisMode::TIME;
static std::string g_plot_title;
static std::vector<Button> toolbar_buttons;


constexpr unsigned long COLOR_BG = 0x000000; // Black background
constexpr unsigned long COLOR_FG = 0xFFFFFF; // White foreground
constexpr unsigned long COLOR_LINE = 0x2CFF05; // Green line color

constexpr double LEFT_MARGIN = 0.08;
constexpr double RIGHT_MARGIN = 0.03;
constexpr double TOP_MARGIN = 0.10;
constexpr double BOTTOM_MARGIN = 0.15;
constexpr int TOOLBAR_HEIGHT = 30;
constexpr int LEGEND_HEIGHT = 16;

static int g_line_thickness = 1;

void render_thick_point(Display* dpy, Drawable drawable, GC gc, int x, int y, int radius) {
    int r = std::max(1, radius);
    XFillArc(dpy, drawable, gc, x - r, y - r, 2 * r, 2 * r, 0, 360 * 64);
}

static std::vector<PlotSample> process_buffer(
    const void* data, std::size_t num_elements, std::size_t element_size_bytes,
    bool is_complex, bool is_float, double xstart, double xdelta
) {
    std::vector<PlotSample> samples;
    if (is_float) {
        const float* fdata = static_cast<const float*>(data);
        for (std::size_t i = 0; i < num_elements; ++i) {
            float re = fdata[is_complex ? 2*i : i];
            float im = is_complex ? fdata[2*i + 1] : 0;
            samples.push_back({xstart + i * xdelta, re, im});
        }
    } else {
        const short* idata = static_cast<const short*>(data);
        for (std::size_t i = 0; i < num_elements; ++i) {
            float re = idata[is_complex ? 2*i : i];
            float im = is_complex ? idata[2*i + 1] : 0;
            samples.push_back({xstart + i * xdelta, re, im});
        }
    }
    return samples;
}

static double compute_value(const PlotSample& s, PlotMode mode) {
    switch (mode) {
        case PlotMode::Real: return s.real;
        case PlotMode::Imag: return s.imag;
        case PlotMode::Magnitude: return std::sqrt(s.real * s.real + s.imag * s.imag);
        case PlotMode::Phase: return std::atan2(s.real, s.imag);
    }
    return s.real;
}

static std::string mode_to_string(PlotMode mode) {
    switch (mode) {
        case PlotMode::Magnitude: return "Magnitude";
        case PlotMode::Real: return "Real";
        case PlotMode::Imag: return "Imag";
        case PlotMode::Phase: return "Phase";
        default: return "Unknown";
    }
}


struct DecimatedTrace {
    std::vector<std::pair<double, double>> y_minmax_per_pixel;
    ZoomRegion cached_region;
    PlotMode cached_mode;
    int cached_width;
    bool dirty = true;

    void update_if_needed(const std::vector<PlotSample>& samples, PlotMode mode,
                          const ZoomRegion& view, int pixel_width, bool use_index) {
        if (!dirty && cached_region.xmin == view.xmin && cached_region.xmax == view.xmax &&
            cached_mode == mode && cached_width == pixel_width) return;

        y_minmax_per_pixel.assign(pixel_width, {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()});

        for (size_t i = 0; i < samples.size(); ++i) {
            double xval = use_index ? (double)i : samples[i].time;
            if (xval < view.xmin || xval > view.xmax) continue;

            int bin = (int)((xval - view.xmin) / (view.xmax - view.xmin) * pixel_width);
            if (bin < 0 || bin >= pixel_width) continue;

            double val = compute_value(samples[i], mode);
            auto& [minv, maxv] = y_minmax_per_pixel[bin];
            if (val < minv) minv = val;
            if (val > maxv) maxv = val;
        }

        cached_mode = mode;
        cached_region = view;
        cached_width = pixel_width;
        dirty = false;
    }
};

static std::vector<DecimatedTrace> g_cached_decimated_traces;


static ZoomRegion autoscale_region(const std::vector<PlotSample>& samples, PlotMode mode, std::optional<std::pair<double, double>> y_range = std::nullopt) {
    double tmin, tmax;

    if (g_xaxis_mode == XAxisMode::INDEX) {
        tmin = 0;
        tmax = samples.size() - 1;
    } else {
        tmin = samples.front().time;
        tmax = samples.back().time;
    }

    double vmin = compute_value(samples.front(), mode);
    double vmax = vmin;
    if (y_range) {
        vmin = y_range->first;
        vmax = y_range->second;
    } else {
        for (const auto& s : samples) {
            double v = compute_value(s, mode);
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
    }

    double vrange = vmax - vmin;
    if (vrange == 0) {
        vmin -= 1.0; vmax += 1.0;
    } else {
        if (!y_range) {
            vmin -= 0.05 * vrange;
            vmax += 0.05 * vrange;
        }
    }

    return {tmin, tmax, vmin, vmax};
}

static void save_pixmap_to_ppm(Display* dpy, Drawable pixmap, int width, int height, const char* filename) {
    XImage* image = XGetImage(dpy, pixmap, 0, 0, width, height, AllPlanes, ZPixmap);
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            unsigned long p = XGetPixel(image, x, y);
            ofs.put((p >> 16) & 0xFF); ofs.put((p >> 8) & 0xFF); ofs.put(p & 0xFF);
        }
    XDestroyImage(image);
    std::cout << "Saved image to " << filename << "\n";
}

static void draw_text(Display* dpy, Drawable win, GC gc, int x, int y, const std::string& str) {
    XDrawString(dpy, win, gc, x, y, str.c_str(), str.length());
}

static void render_axes(Display* dpy, Drawable win, GC gc,
    const ZoomRegion& r, int px, int py, int pw, int ph,
    int num_samples, bool use_index) {
    const int ticks = 5;

    for (int i = 0; i <= ticks; ++i) {
        double val = r.xmin + i * (r.xmax - r.xmin) / ticks;

        int x = px + i * pw / ticks;
        XDrawLine(dpy, win, gc, x, py + ph, x, py + ph + 5);

        char label[32];
        if (use_index)
            snprintf(label, sizeof(label), "%d", (int)val);
        else
            snprintf(label, sizeof(label), "%.2f", val);

        draw_text(dpy, win, gc, x - 10, py + ph + 20, label);
    }

    for (int i = 0; i <= ticks; ++i) {
        double v = r.ymin + i * (r.ymax - r.ymin) / ticks;
        int y = py + (int)((r.ymax - v) / (r.ymax - r.ymin) * ph);

        XDrawLine(dpy, win, gc, px - 5, y, px, y);

        std::ostringstream ss;
        ss.precision(2);
        ss << std::fixed << v;
        std::string label = ss.str();

        int text_width = XTextWidth(XQueryFont(dpy, XGContextFromGC(gc)), label.c_str(), label.length());
        draw_text(dpy, win, gc, px - 8 - text_width, y + 4, label);
    }
}

void prepare_decimated_traces(const std::vector<std::vector<PlotSample>>& traces,
                              PlotMode mode, const ZoomRegion& view,
                              int pixel_width, bool use_index) {
    if (g_cached_decimated_traces.size() != traces.size())
        g_cached_decimated_traces.resize(traces.size());

    for (std::size_t i = 0; i < traces.size(); ++i) {
        g_cached_decimated_traces[i].update_if_needed(traces[i], mode, view, pixel_width, use_index);
    }
}

// Part 3/3: Optimized render_pixmap()

void render_pixmap(Display* dpy, Pixmap pixmap, GC gc, int w, int h,
                   const ZoomRegion& r, const std::string& title) {
    XSetForeground(dpy, gc, COLOR_BG);
    XFillRectangle(dpy, pixmap, gc, 0, 0, w, h);

    int plot_width = w * 0.9;
    int plot_height = h * 0.72;
    int px = (w - plot_width) / 2;
    int py = (h - plot_height - TOOLBAR_HEIGHT - LEGEND_HEIGHT - 10) / 2;
    int pw = plot_width;
    int ph = plot_height;

    bool use_index = g_xaxis_mode == XAxisMode::INDEX;
    prepare_decimated_traces(g_traces, g_mode, r, pw, use_index);

    std::vector<unsigned long> trace_colors = {
        0x2CFF05, 0xFF0000, 0x00C0FF, 0xFF00FF, 0xFFFF00,
        0xFFA500, 0x00FFFF, 0xFFFFFF, 0xFF69B4, 0xB22222
    };

    XGCValues old_vals;
    XGetGCValues(dpy, gc, GCLineWidth | GCLineStyle | GCCapStyle | GCJoinStyle, &old_vals);
    XSetLineAttributes(dpy, gc, g_line_thickness, LineSolid, CapButt, JoinMiter);

    // Draw all traces
    for (std::size_t t = 0; t < g_traces.size(); ++t) {
        if (!g_trace_visibility[t]) continue;
        if (t >= g_cached_decimated_traces.size()) continue;
        const auto& bins = g_cached_decimated_traces[t].y_minmax_per_pixel;
        if (bins.empty()) continue;

        XSetForeground(dpy, gc, trace_colors[t % trace_colors.size()]);
        int prev_x = -1, prev_y = -1;

        for (int xpix = 0; xpix < (int)bins.size(); ++xpix) {
            const auto& [minv, maxv] = bins[xpix];
            if (minv > maxv) continue;

            int screen_x = px + xpix;
            int y1 = py + (int)((r.ymax - minv) / (r.ymax - r.ymin) * ph);
            int y2 = py + (int)((r.ymax - maxv) / (r.ymax - r.ymin) * ph);
            y1 = std::clamp(y1, py, py + ph - 1);
            y2 = std::clamp(y2, py, py + ph - 1);

            if (g_trace_styles[t] == PlotRenderStyle::LINES) {
                XDrawLine(dpy, pixmap, gc, screen_x, y1, screen_x, y2);
                if (prev_x != -1 && std::abs(screen_x - prev_x) >= 1) {
                    XDrawLine(dpy, pixmap, gc, prev_x, prev_y, screen_x, y1);
                }
                prev_x = screen_x;
                prev_y = y2;
            } else {
                render_thick_point(dpy, pixmap, gc, screen_x, y1, g_line_thickness);
                if (minv != maxv)
                    render_thick_point(dpy, pixmap, gc, screen_x, y2, g_line_thickness);
            }
        }
    }

    // Restore GC before drawing axes and borders
    XSetForeground(dpy, gc, COLOR_FG);
    XSetLineAttributes(dpy, gc, 1, LineSolid, CapButt, JoinMiter);

    // Draw white border box around plot
    XDrawRectangle(dpy, pixmap, gc, px, py, pw, ph);

    // Axes
    render_axes(dpy, pixmap, gc, r, px, py, pw, ph, g_traces[0].size(), use_index);

    // Plot title
    draw_text(dpy, pixmap, gc, w / 2 - (title.length() * 3), py - 10, title);


    g_legend_boxes.clear();  // Reset each frame

    int legend_start_x = px + 10;
    int legend_y = h - TOOLBAR_HEIGHT - LEGEND_HEIGHT - 30;
    int box_height = 14;
    int box_width = 12;
    int spacing_x = 100;

    for (std::size_t t = 0; t < g_traces.size(); ++t) {
        int lx = legend_start_x + t * spacing_x;

        // Dim if hidden
        if (!g_trace_visibility[t]) {
            XSetForeground(dpy, gc, 0x404040);  // dark gray
        } else {
            XSetForeground(dpy, gc, trace_colors[t % trace_colors.size()]);
        }

        XFillRectangle(dpy, pixmap, gc, lx, legend_y, box_width, box_height);

        XSetForeground(dpy, gc, COLOR_FG);
        std::string label = "Trace " + std::to_string(t);
        draw_text(dpy, pixmap, gc, lx + box_width + 6, legend_y + 10, label);

        g_legend_boxes.push_back({lx, legend_y, spacing_x, box_height});
    }

    // Toolbar
    int tool_y = h - TOOLBAR_HEIGHT, xpos = 20;
    for (auto& b : toolbar_buttons) {
        XDrawRectangle(dpy, pixmap, gc, xpos, tool_y, 100, TOOLBAR_HEIGHT - 10);
        draw_text(dpy, pixmap, gc, xpos + 10, tool_y + 18, b.label);
        b.x = xpos; b.width = 100;
        xpos += 120;
    }
}




// Continuation of the plot_buffer and event loop with Toggle Style button support

void plot_buffer(const void* data, std::size_t num_elements, std::size_t elem_bytes,
    bool is_complex, bool is_float, double xstart, double xdelta,
    const std::string& plot_title, int line_thickness,
    std::optional<std::pair<double, double>> y_range,
    std::size_t num_traces) {

    g_plot_title = plot_title;
    g_line_thickness = line_thickness;
    g_traces.clear();

    g_trace_styles.clear();
    g_trace_styles.resize(num_traces, PlotRenderStyle::LINES);
    
    g_trace_visibility.clear();
    g_trace_visibility.resize(num_traces, true);

    if (num_traces == 0) return;
    std::size_t trace_len = num_elements / num_traces;

    const char* base = static_cast<const char*>(data);
    trace_len = std::min(trace_len, num_elements);  // Clamp to avoid overrun
    
    int elem_bytes_complex = is_complex ? 2 * elem_bytes : elem_bytes;
    for (std::size_t t = 0; t < num_traces; ++t) {
        const void* trace_ptr = base + t * trace_len * elem_bytes_complex;
        auto samples = process_buffer(trace_ptr, trace_len, elem_bytes, is_complex, is_float, xstart, xdelta);
        if (!samples.empty()) {
            g_traces.push_back(std::move(samples));
        }
    }
    
    if (g_traces.empty()) {
        std::cerr << "No valid traces found.\n";
        return;
    }

    Display* dpy = XOpenDisplay(nullptr);
    if (!dpy) {
        std::cerr << "X11 not available.\n";
        return;
    }

    int screen = DefaultScreen(dpy);
    int w = 1000, h = 600;
    Window win = XCreateSimpleWindow(dpy, RootWindow(dpy, screen), 10, 10, w, h, 1,
                                     BlackPixel(dpy, screen), WhitePixel(dpy, screen));
    XSelectInput(dpy, win, ExposureMask | ButtonPressMask | ButtonReleaseMask |
                        PointerMotionMask | StructureNotifyMask);
    XMapWindow(dpy, win);
    GC gc = XCreateGC(dpy, win, 0, nullptr);

    Pixmap pixmap = XCreatePixmap(dpy, win, w, h, DefaultDepth(dpy, screen));
    bool pixmap_dirty = true;

    
    
    ZoomRegion view = autoscale_region(g_traces[0], g_mode, y_range);
    std::cout << "Plotting " << g_traces.size() << " traces with mode: " << mode_to_string(g_mode) << "\n";

    int bx0 = -1, by0 = -1, bx1 = -1, by1 = -1;
    bool dragging = false;
    std::string readout;

    toolbar_buttons = {
        {"Save PNG", 0, 100},
        {"Cycle Mode", 0, 100},
        {"Reset Zoom", 0, 100},
        {"Cycle X-Axis", 0, 120}
    };


    XEvent e;
    while (true) {
        if (XPending(dpy)) {
            XNextEvent(dpy, &e);
            switch (e.type) {
                case ConfigureNotify:
                    if (e.xconfigure.width != w || e.xconfigure.height != h) {
                        w = e.xconfigure.width;
                        h = e.xconfigure.height;
                        XFreePixmap(dpy, pixmap);
                        pixmap = XCreatePixmap(dpy, win, w, h, DefaultDepth(dpy, screen));
                        pixmap_dirty = true;
                    }
                    break;

                case MotionNotify: {
                    int mx = e.xmotion.x;
                    int my = e.xmotion.y;

                    int plot_width = w * 0.9;
                    int plot_height = h * 0.72;
                    int px = (w - plot_width) / 2;
                    int py = (h - plot_height - TOOLBAR_HEIGHT - LEGEND_HEIGHT - 10) / 2;
                    int pw = plot_width;
                    int ph = plot_height;

                    // Clamp mouse position to plot box
                    int rel_x = std::clamp(mx - px, 0, pw - 1);
                    int rel_y = std::clamp(my - py, 0, ph - 1);

                    // Use exact screen-to-data conversion
                    double t = view.xmin + (rel_x / (double)pw) * (view.xmax - view.xmin);
                    double v = view.ymax - (rel_y / (double)ph) * (view.ymax - view.ymin);

                    std::ostringstream ss;
                    ss.precision(2);
                    if (g_xaxis_mode == XAxisMode::TIME)
                        ss << "Time: " << t;
                    else {
                        int display_idx = std::clamp(static_cast<int>(std::round(t)), 0, static_cast<int>(g_traces[0].size()) - 1);
                        ss << "Index: " << display_idx;
                    }
                    ss << "  Value: " << v;
                    readout = ss.str();

                    if (dragging) {
                        bx1 = mx;
                        by1 = my;
                    }

                    XCopyArea(dpy, pixmap, win, gc, 0, 0, w, h, 0, 0);

                    int label_y = h - TOOLBAR_HEIGHT - 10;
                    XSetForeground(dpy, gc, COLOR_FG);
                    draw_text(dpy, win, gc, 40, label_y, "Mode: " + mode_to_string(g_mode) + "  " + readout);

                    if (bx0 >= 0 && bx1 >= 0) {
                        int zx = std::min(bx0, bx1), zy = std::min(by0, by1);
                        int zw = std::abs(bx1 - bx0), zh = std::abs(by1 - by0);
                        XSetForeground(dpy, gc, COLOR_FG);
                        XDrawRectangle(dpy, win, gc, zx, zy, zw, zh);
                    }

                    XFlush(dpy);
                    break;
                }

                case ButtonPress:
                    if (e.xbutton.y >= h - TOOLBAR_HEIGHT - 10) {
                        for (auto& b : toolbar_buttons) {
                            if (e.xbutton.x >= b.x && e.xbutton.x <= b.x + b.width) {
                                if (b.label == "Save PNG") {
                                    save_pixmap_to_ppm(dpy, pixmap, w, h, "plot_out.ppm");
                                } else {
                                    auto y_range_tmp = y_range;
                                    if (b.label == "Cycle Mode") {
                                        g_mode = static_cast<PlotMode>((static_cast<int>(g_mode) + 1) % 4);
                                        if (g_mode == PlotMode::Phase) y_range_tmp = std::make_pair(-M_PI, M_PI);
                                    } else if (b.label == "Reset Zoom") {
                                        // Placeholder until we figure out a way to not reset zoom on all toolbar changes
                                        // Likely means isolating our auto scale to x and y axis, 
                                    } else if (b.label == "Cycle X-Axis") {
                                        g_xaxis_mode = (g_xaxis_mode == XAxisMode::TIME) ? XAxisMode::INDEX : XAxisMode::TIME;
                                    }

                                    view = autoscale_region(g_traces[0], g_mode, y_range);  // <-- Recalculate view bounds
                                    while (!g_zoom_stack.empty()) g_zoom_stack.pop();
                                    g_zoom_depth = 0;
                                    pixmap_dirty = true;
                                }
                            }
                        }
                    } else if (e.xbutton.y >= h - TOOLBAR_HEIGHT - LEGEND_HEIGHT - 30) {
                        for (std::size_t t = 0; t < g_legend_boxes.size(); ++t) {
                            const auto& box = g_legend_boxes[t];
                            if (e.xbutton.y >= box.y && e.xbutton.y <= box.y + box.height &&
                                e.xbutton.x >= box.x && e.xbutton.x <= box.x + box.width) {
    
                                if (e.xbutton.button == Button1) {
                                    g_trace_visibility[t] = !g_trace_visibility[t];
                                } else if (e.xbutton.button == Button3) {
                                    g_trace_styles[t] = (g_trace_styles[t] == PlotRenderStyle::LINES)
                                        ? PlotRenderStyle::DOTS
                                        : PlotRenderStyle::LINES;
                                }
                                pixmap_dirty = true;
                            }
                        }
                    } else {
                        if (e.xbutton.button == Button1) {
                            dragging = true;
                            bx0 = bx1 = e.xbutton.x;
                            by0 = by1 = e.xbutton.y;
                        } else if (e.xbutton.button == Button2) {
                            g_mode = static_cast<PlotMode>((static_cast<int>(g_mode) + 1) % 4);
                            view = autoscale_region(g_traces[0], g_mode);
                            while (!g_zoom_stack.empty()) g_zoom_stack.pop();
                            g_zoom_depth = 0;
                            pixmap_dirty = true;
                        } else if (e.xbutton.button == Button3) {
                            if (!g_zoom_stack.empty()) {
                                view = g_zoom_stack.top();
                                g_zoom_stack.pop();
                                g_zoom_depth--;
                                pixmap_dirty = true;
                            }
                        }
                    }
                    break;

                case ButtonRelease:
                    if (dragging) {
                        dragging = false;
                        int x1 = std::min(bx0, bx1), x2 = std::max(bx0, bx1);
                        int y1 = std::min(by0, by1), y2 = std::max(by0, by1);
                        int px = w * LEFT_MARGIN, py = h * TOP_MARGIN;
                        int pw = w * (1 - LEFT_MARGIN - RIGHT_MARGIN);
                        int ph = h * (1 - TOP_MARGIN - BOTTOM_MARGIN);

                        
                        if (g_zoom_depth < MAX_ZOOM_HISTORY) {
                            if (x2 - x1 > 10 && y2 - y1 > 10) {
                                g_zoom_stack.push(view);
                                g_zoom_depth++;
                                double nxmin = view.xmin + (x1 - px) / (double)pw * (view.xmax - view.xmin);
                                double nxmax = view.xmin + (x2 - px) / (double)pw * (view.xmax - view.xmin);
                                double nymin = view.ymax - (y2 - py) / (double)ph * (view.ymax - view.ymin);
                                double nymax = view.ymax - (y1 - py) / (double)ph * (view.ymax - view.ymin);
                                view = {nxmin, nxmax, nymin, nymax};
                                pixmap_dirty = true;
                            }
                        }

                        bx0 = bx1 = by0 = by1 = -1;
                    }
                    break;
            }
        } else {
            if (pixmap_dirty) {
                render_pixmap(dpy, pixmap, gc, w, h, view, g_plot_title);
                pixmap_dirty = false;
            }

            XCopyArea(dpy, pixmap, win, gc, 0, 0, w, h, 0, 0);

            int label_y = h - TOOLBAR_HEIGHT - 10;
            XSetForeground(dpy, gc, COLOR_FG);
            draw_text(dpy, win, gc, 40, label_y, "Mode: " + mode_to_string(g_mode) + "  " + readout);

            if (bx0 >= 0 && bx1 >= 0) {
                int zx = std::min(bx0, bx1), zy = std::min(by0, by1);
                int zw = std::abs(bx1 - bx0), zh = std::abs(by1 - by0);
                XSetForeground(dpy, gc, COLOR_FG);
                XDrawRectangle(dpy, win, gc, zx, zy, zw, zh);
            }

            XFlush(dpy);
            usleep(10000);
        }
    }

    XFreePixmap(dpy, pixmap);
    XFreeGC(dpy, gc);
    XDestroyWindow(dpy, win);
    XCloseDisplay(dpy);
}

}  // namespace xplot