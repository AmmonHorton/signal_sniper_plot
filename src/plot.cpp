#include "plot.hpp"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <stack>
#include <sstream>
#include <cstring>

namespace xplot {

enum class PlotMode { Magnitude, Real, Imag, Phase };

struct PlotSample {
    double time;
    double real;
    double imag;
};

struct ZoomRegion {
    double xmin, xmax;
    double ymin, ymax;
};

static std::vector<PlotSample> g_samples;
static std::stack<ZoomRegion> g_zoom_stack;
static PlotMode g_mode = PlotMode::Magnitude;

// ---- Configurable colors ----
static unsigned long COLOR_BG, COLOR_FG, COLOR_BOX, COLOR_TEXT;

static std::vector<PlotSample> process_buffer(
    const void* data,
    std::size_t num_elements,
    std::size_t element_size_bytes,
    bool is_complex,
    bool is_float,
    double xstart,
    double xdelta
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
        case PlotMode::Phase: return std::atan2(s.imag, s.real);
        default: return s.real;
    }
}

static ZoomRegion autoscale_region(const std::vector<PlotSample>& samples, PlotMode mode) {
    double tmin = samples.front().time;
    double tmax = samples.back().time;
    double vmin = compute_value(samples.front(), mode);
    double vmax = vmin;

    for (const auto& s : samples) {
        double v = compute_value(s, mode);
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }

    // Pad a little
    double vrange = vmax - vmin;
    if (vrange == 0) { vmin -= 1.0; vmax += 1.0; }
    else { vmin -= 0.05 * vrange; vmax += 0.05 * vrange; }

    return {tmin, tmax, vmin, vmax};
}

static void draw_text(Display* dpy, Window win, GC gc, int x, int y, const std::string& str) {
    XDrawString(dpy, win, gc, x, y, str.c_str(), str.length());
}

static void render_axes(Display* dpy, Window win, GC gc, const ZoomRegion& r, int width, int height) {
    // Tick marks on X
    int tick_count = 10;
    for (int i = 0; i <= tick_count; ++i) {
        double t = r.xmin + i * (r.xmax - r.xmin) / tick_count;
        int x = static_cast<int>((t - r.xmin) / (r.xmax - r.xmin) * width);
        XDrawLine(dpy, win, gc, x, height - 20, x, height - 10);

        std::ostringstream label;
        label.precision(2); label << std::fixed << t;
        draw_text(dpy, win, gc, x - 15, height - 5, label.str());
    }

    // Tick marks on Y
    for (int i = 0; i <= tick_count; ++i) {
        double v = r.ymin + i * (r.ymax - r.ymin) / tick_count;
        int y = static_cast<int>((r.ymax - v) / (r.ymax - r.ymin) * height);
        XDrawLine(dpy, win, gc, 0, y, 10, y);

        std::ostringstream label;
        label.precision(2); label << std::fixed << v;
        draw_text(dpy, win, gc, 12, y + 4, label.str());
    }
}

static void render_plot(Display* dpy, Window win, GC gc, int width, int height, const ZoomRegion& region, const std::string& status_text = "") {
    XClearWindow(dpy, win);

    std::vector<XPoint> points;
    for (const auto& s : g_samples) {
        if (s.time < region.xmin || s.time > region.xmax) continue;
        double val = compute_value(s, g_mode);
        int x = static_cast<int>((s.time - region.xmin) / (region.xmax - region.xmin) * width);
        int y = static_cast<int>((region.ymax - val) / (region.ymax - region.ymin) * height);
        points.push_back({static_cast<short>(x), static_cast<short>(y)});
    }

    if (points.size() > 1) {
        XSetForeground(dpy, gc, COLOR_FG);
        XDrawLines(dpy, win, gc, points.data(), points.size(), CoordModeOrigin);
    }

    XSetForeground(dpy, gc, COLOR_TEXT);
    render_axes(dpy, win, gc, region, width, height);
    draw_text(dpy, win, gc, 20, 20, status_text);
}

void plot_buffer(
    const void* data,
    std::size_t num_elements,
    std::size_t element_size_bytes,
    bool is_complex,
    bool is_float,
    double xstart_j1950_sec,
    double xdelta_sec
) {
    g_samples = process_buffer(data, num_elements, element_size_bytes, is_complex, is_float, xstart_j1950_sec, xdelta_sec);
    if (g_samples.empty()) return;

    Display* dpy = XOpenDisplay(NULL);
    if (!dpy) {
        std::cerr << "Unable to open X display\n";
        return;
    }

    int screen = DefaultScreen(dpy);
    int width = 1000, height = 600;

    Window win = XCreateSimpleWindow(dpy, RootWindow(dpy, screen), 100, 100, width, height, 1,
                                     BlackPixel(dpy, screen), WhitePixel(dpy, screen));
    XSelectInput(dpy, win, ExposureMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask);
    XMapWindow(dpy, win);
    GC gc = XCreateGC(dpy, win, 0, nullptr);

    // Init colors
    COLOR_BG = WhitePixel(dpy, screen);
    COLOR_FG = BlackPixel(dpy, screen);
    COLOR_BOX = 0xFF0000;  // Red
    COLOR_TEXT = 0x003366; // Navy

    ZoomRegion current = autoscale_region(g_samples, g_mode);

    bool dragging = false;
    int box_x0 = 0, box_y0 = 0, box_x1 = 0, box_y1 = 0;
    std::string mouse_readout = "";

    while (true) {
        XEvent e;
        XNextEvent(dpy, &e);

        if (e.type == Expose) {
            render_plot(dpy, win, gc, width, height, current, mouse_readout);
        } else if (e.type == MotionNotify) {
            double t = current.xmin + (double)e.xmotion.x / width * (current.xmax - current.xmin);
            double v = current.ymax - (double)e.xmotion.y / height * (current.ymax - current.ymin);

            std::ostringstream msg;
            msg.precision(2);
            msg << "Time: " << std::fixed << t << ", Value: " << v;
            mouse_readout = msg.str();
            render_plot(dpy, win, gc, width, height, current, mouse_readout);
        } else if (e.type == ButtonPress) {
            if (e.xbutton.button == Button1) {
                dragging = true;
                box_x0 = e.xbutton.x;
                box_y0 = e.xbutton.y;
            } else if (e.xbutton.button == Button2) {
                g_mode = static_cast<PlotMode>((static_cast<int>(g_mode) + 1) % 4);
                current = autoscale_region(g_samples, g_mode);
                g_zoom_stack = {};
                render_plot(dpy, win, gc, width, height, current, mouse_readout);
            } else if (e.xbutton.button == Button3) {
                if (!g_zoom_stack.empty()) {
                    current = g_zoom_stack.top();
                    g_zoom_stack.pop();
                    render_plot(dpy, win, gc, width, height, current, mouse_readout);
                }
            }
        } else if (e.type == ButtonRelease && dragging && e.xbutton.button == Button1) {
            dragging = false;
            box_x1 = e.xbutton.x;
            box_y1 = e.xbutton.y;

            int x_start = std::min(box_x0, box_x1);
            int x_end   = std::max(box_x0, box_x1);
            int y_start = std::min(box_y0, box_y1);
            int y_end   = std::max(box_y0, box_y1);

            if ((x_end - x_start > 10) && (y_end - y_start > 10)) {
                g_zoom_stack.push(current);
                double new_xmin = current.xmin + x_start / double(width) * (current.xmax - current.xmin);
                double new_xmax = current.xmin + x_end   / double(width) * (current.xmax - current.xmin);
                double new_ymax = current.ymin + y_start / double(height) * (current.ymax - current.ymin);
                double new_ymin = current.ymin + y_end   / double(height) * (current.ymax - current.ymin);
                current = {new_xmin, new_xmax, new_ymin, new_ymax};
                render_plot(dpy, win, gc, width, height, current, mouse_readout);
            }
        }
    }

    XCloseDisplay(dpy);
}

}  // namespace xplot
