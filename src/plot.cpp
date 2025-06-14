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

enum class XAxisMode { TIME, INDEX };

static XAxisMode g_xaxis_mode = XAxisMode::TIME;
static std::string g_plot_title;

static std::vector<PlotSample> g_samples;
static std::stack<ZoomRegion> g_zoom_stack;
static PlotMode g_mode = PlotMode::Magnitude;

static unsigned long COLOR_BG, COLOR_FG, COLOR_BOX, COLOR_TEXT;

constexpr double LEFT_MARGIN = 0.08;
constexpr double RIGHT_MARGIN = 0.03;
constexpr double TOP_MARGIN = 0.10;
constexpr double BOTTOM_MARGIN = 0.15;
constexpr int TOOLBAR_HEIGHT = 30;

struct Button { std::string label; int x, width; };
std::vector<Button> toolbar_buttons;

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
        case PlotMode::Phase: return std::atan2(s.imag, s.real);
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
    double vrange = vmax - vmin;
    if (vrange == 0) { vmin -= 1.0; vmax += 1.0; }
    else { vmin -= 0.05 * vrange; vmax += 0.05 * vrange; }
    return {tmin, tmax, vmin, vmax};
}

static void save_window_to_ppm(Display* dpy, Window win, int width, int height, const char* filename) {
    XImage* image = XGetImage(dpy, win, 0, 0, width, height, AllPlanes, ZPixmap);
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

// X-axis ticks and labels
for (int i = 0; i <= ticks; ++i) {
double val = use_index
? i * (num_samples - 1) / (double)ticks
: r.xmin + i * (r.xmax - r.xmin) / ticks;

int x = px + i * pw / ticks;
XDrawLine(dpy, win, gc, x, py + ph, x, py + ph + 5);

char label[32];
if (use_index)
snprintf(label, sizeof(label), "%d", (int)val);
else
snprintf(label, sizeof(label), "%.2f", val);

draw_text(dpy, win, gc, x - 10, py + ph + 20, label);
}

// Y-axis ticks and right-aligned labels
for (int i = 0; i <= ticks; ++i) {
double v = r.ymin + i * (r.ymax - r.ymin) / ticks;
int y = py + (int)((r.ymax - v) / (r.ymax - r.ymin) * ph);

XDrawLine(dpy, win, gc, px - 5, y, px, y);

std::ostringstream ss;
ss.precision(2);
ss << std::fixed << v;
std::string label = ss.str();

// Measure label width to right-align it
int text_width = XTextWidth(XQueryFont(dpy, XGContextFromGC(gc)), label.c_str(), label.length());
draw_text(dpy, win, gc, px - 8 - text_width, y + 4, label);  // right-aligned
}
}

static void render_plot(Display* dpy, Window win, GC gc, int w, int h,
    const ZoomRegion& r, const std::string& readout,
    int bx0 = -1, int by0 = -1, int bx1 = -1, int by1 = -1)
{
    XClearWindow(dpy, win);

    // Centered plot box (90% width, 85% height)
    int plot_width = w * 0.9;
    int plot_height = h * 0.75;
    int px = (w - plot_width) / 2;
    int py = (h - plot_height - TOOLBAR_HEIGHT - 30) / 2;  // 30 px buffer for title
    int pw = plot_width;
    int ph = plot_height;

    // Draw signal inside clipping box
    std::vector<XPoint> pts;
    for (auto& s : g_samples) {
        if (s.time < r.xmin || s.time > r.xmax) continue;
        double val = compute_value(s, g_mode);
        int x = px + (int)((s.time - r.xmin) / (r.xmax - r.xmin) * pw);
        int y = py + (int)((r.ymax - val) / (r.ymax - r.ymin) * ph);
        pts.push_back({(short)x, (short)y});
    }

    if (pts.size() > 1) {
        XRectangle clip;
        clip.x = px;
        clip.y = py;
        clip.width = pw;
        clip.height = ph;
        XSetClipRectangles(dpy, gc, 0, 0, &clip, 1, Unsorted);

        XSetForeground(dpy, gc, COLOR_FG);
        XDrawLines(dpy, win, gc, pts.data(), pts.size(), CoordModeOrigin);

        XSetClipMask(dpy, gc, None);
    }

    // Bounding box
    XSetForeground(dpy, gc, COLOR_FG);
    XDrawRectangle(dpy, win, gc, px, py, pw, ph);

    // Axis ticks and labels
    render_axes(dpy, win, gc, r, px, py, pw, ph, g_samples.size(), g_xaxis_mode == XAxisMode::INDEX);

    // Title
    draw_text(dpy, win, gc, w / 2 - (g_plot_title.length() * 3), py - 10, g_plot_title);

    // Mode and readout (below X-axis)
    draw_text(dpy, win, gc, px, py + ph + 40, "Mode: " + mode_to_string(g_mode));
    draw_text(dpy, win, gc, px + 200, py + ph + 40, readout);

    // Toolbar
    int tool_y = h - TOOLBAR_HEIGHT, xpos = 20;
    for (auto& b : toolbar_buttons) {
        XDrawRectangle(dpy, win, gc, xpos, tool_y, 100, TOOLBAR_HEIGHT - 10);
        draw_text(dpy, win, gc, xpos + 10, tool_y + 18, b.label);
        b.x = xpos; b.width = 100;
        xpos += 120;
    }

    // Zoom box
    if (bx0 >= 0 && bx1 >= 0) {
        int zx = std::min(bx0, bx1), zy = std::min(by0, by1);
        int zw = std::abs(bx1 - bx0), zh = std::abs(by1 - by0);
        XSetForeground(dpy, gc, COLOR_BOX);
        XDrawRectangle(dpy, win, gc, zx, zy, zw, zh);
    }

    XFlush(dpy);
}

void plot_buffer(const void* data, std::size_t num_elements, std::size_t elem_bytes,
bool is_complex, bool is_float, double xstart, double xdelta,
const std::string& plot_title) {
g_plot_title = plot_title;
g_samples = process_buffer(data, num_elements, elem_bytes, is_complex, is_float, xstart, xdelta);
if (g_samples.empty()) return;

Display* dpy = XOpenDisplay(nullptr);
if (!dpy) { std::cerr << "X11 not available.\n"; return; }

int screen = DefaultScreen(dpy);
int w = 1000, h = 600;
Window win = XCreateSimpleWindow(dpy, RootWindow(dpy, screen), 10, 10, w, h, 1,
                 BlackPixel(dpy, screen), WhitePixel(dpy, screen));
XSelectInput(dpy, win, ExposureMask | ButtonPressMask | ButtonReleaseMask |
        PointerMotionMask | StructureNotifyMask);
XMapWindow(dpy, win);
GC gc = XCreateGC(dpy, win, 0, nullptr);

COLOR_BG = WhitePixel(dpy, screen);
COLOR_FG = BlackPixel(dpy, screen);
COLOR_BOX = 0xFF0000;
COLOR_TEXT = 0x003366;

ZoomRegion view = autoscale_region(g_samples, g_mode);
int bx0 = -1, by0 = -1, bx1 = -1, by1 = -1;
bool dragging = false;
std::string readout;

toolbar_buttons = {
{"Save PNG", 0, 100},
{"Cycle Mode", 0, 100},
{"Reset Zoom", 0, 100},
{"Cycle X-Axis", 0, 120}
};
bool needs_redraw = true;
XEvent e;

while (true) {
    if (XPending(dpy)) {
        XNextEvent(dpy, &e);
        switch (e.type) {
            case ConfigureNotify:
                w = e.xconfigure.width;
                h = e.xconfigure.height;
                needs_redraw = true;
                break;

            case Expose:
                needs_redraw = true;
                break;

            case MotionNotify: {
                int mx = e.xmotion.x;
                int my = e.xmotion.y;

                double t = (g_xaxis_mode == XAxisMode::TIME)
                    ? view.xmin + (mx - w * LEFT_MARGIN) / (w * (1 - LEFT_MARGIN - RIGHT_MARGIN)) * (view.xmax - view.xmin)
                    : (mx - w * LEFT_MARGIN) / (w * (1 - LEFT_MARGIN - RIGHT_MARGIN)) * (g_samples.size() - 1);

                double v = view.ymax - (my - h * TOP_MARGIN) / (h * (1 - TOP_MARGIN - BOTTOM_MARGIN)) * (view.ymax - view.ymin);

                std::ostringstream ss;
                ss.precision(2);
                if (g_xaxis_mode == XAxisMode::TIME)
                    ss << "Time: " << t;
                else
                    ss << "Index: " << (int)t;
                ss << "  Value: " << v;

                if (readout != ss.str()) {
                    readout = ss.str();
                    needs_redraw = true;
                }

                if (dragging) {
                    bx1 = mx;
                    by1 = my;
                    needs_redraw = true;
                }
                break;
            }

            case ButtonPress: {
                if (e.xbutton.y >= h - TOOLBAR_HEIGHT) {
                    for (auto& b : toolbar_buttons) {
                        if (e.xbutton.x >= b.x && e.xbutton.x <= b.x + b.width) {
                            if (b.label == "Save PNG")
                                save_window_to_ppm(dpy, win, w, h, "plot_out.ppm");
                            else if (b.label == "Cycle Mode") {
                                g_mode = static_cast<PlotMode>((static_cast<int>(g_mode) + 1) % 4);
                                view = autoscale_region(g_samples, g_mode);
                                while (!g_zoom_stack.empty()) g_zoom_stack.pop();
                                needs_redraw = true;
                            }
                            else if (b.label == "Reset Zoom") {
                                view = autoscale_region(g_samples, g_mode);
                                while (!g_zoom_stack.empty()) g_zoom_stack.pop();
                                needs_redraw = true;
                            }
                            else if (b.label == "Cycle X-Axis") {
                                g_xaxis_mode = (g_xaxis_mode == XAxisMode::TIME)
                                    ? XAxisMode::INDEX
                                    : XAxisMode::TIME;
                                needs_redraw = true;
                            }
                        }
                    }
                } else if (e.xbutton.button == Button1) {
                    dragging = true;
                    bx0 = bx1 = e.xbutton.x;
                    by0 = by1 = e.xbutton.y;
                } else if (e.xbutton.button == Button2) {
                    g_mode = static_cast<PlotMode>((static_cast<int>(g_mode) + 1) % 4);
                    view = autoscale_region(g_samples, g_mode);
                    while (!g_zoom_stack.empty()) g_zoom_stack.pop();
                    needs_redraw = true;
                } else if (e.xbutton.button == Button3) {
                    if (!g_zoom_stack.empty()) {
                        view = g_zoom_stack.top();
                        g_zoom_stack.pop();
                        needs_redraw = true;
                    }
                }
                break;
            }

            case ButtonRelease:
                if (dragging) {
                    dragging = false;
                    int x1 = std::min(bx0, bx1), x2 = std::max(bx0, bx1);
                    int y1 = std::min(by0, by1), y2 = std::max(by0, by1);
                    int px = w * LEFT_MARGIN, py = h * TOP_MARGIN;
                    int pw = w * (1 - LEFT_MARGIN - RIGHT_MARGIN);
                    int ph = h * (1 - TOP_MARGIN - BOTTOM_MARGIN);

                    if (x2 - x1 > 10 && y2 - y1 > 10) {
                        g_zoom_stack.push(view);
                        double nxmin = view.xmin + (x1 - px) / (double)pw * (view.xmax - view.xmin);
                        double nxmax = view.xmin + (x2 - px) / (double)pw * (view.xmax - view.xmin);
                        double nymax = view.ymin + (y1 - py) / (double)ph * (view.ymax - view.ymin);
                        double nymin = view.ymin + (y2 - py) / (double)ph * (view.ymax - view.ymin);
                        view = {nxmin, nxmax, nymin, nymax};
                        needs_redraw = true;
                    }

                    bx0 = bx1 = by0 = by1 = -1;
                }
                break;
        }
    } else {
        // Idle: render once if needed
        if (needs_redraw) {
            render_plot(dpy, win, gc, w, h, view, readout, bx0, by0, bx1, by1);
            needs_redraw = false;
        }
        usleep(500);  // Small pause to avoid CPU spin when idle
    }
}

}

}  // namespace xplot
