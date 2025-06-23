// plot_session.cpp
#include "plot_session.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <unistd.h>

namespace xplot {

constexpr unsigned long COLOR_BG = 0x000000;
constexpr unsigned long COLOR_FG = 0xFFFFFF;

constexpr double LEFT_MARGIN = 0.08;
constexpr double RIGHT_MARGIN = 0.03;
constexpr double TOP_MARGIN = 0.05;
constexpr double BOTTOM_MARGIN = 0.05;
constexpr int READOUT_POS_FROM_RIGHT = 300;
constexpr int TOOLBAR_HEIGHT = 30;
constexpr int VERTICAL_PADDING = 15;
constexpr int LEGEND_HEIGHT = 20;
constexpr int MAX_ZOOM_HISTORY = 5;

static std::vector<unsigned long> trace_colors = {
    0x2CFF05, 0xFF0000, 0x00C0FF, 0xFF00FF, 0xFFFF00,
    0xFFA500, 0x00FFFF, 0xFFFFFF, 0xFF69B4, 0xB22222
};

std::string mode_to_string(PlotMode mode) {
    switch (mode) {
        case PlotMode::Magnitude: return "Magnitude";
        case PlotMode::Real: return "Real";
        case PlotMode::Imag: return "Imag";
        case PlotMode::Phase: return "Phase";
        case PlotMode::IQ: return "I/Q";
        default: return "Unknown";
    }
}

double compute_value(const PlotSample& s, PlotMode mode) {
    switch (mode) {
        case PlotMode::Real: return s.real;
        case PlotMode::Imag: return s.imag;
        case PlotMode::Magnitude: return std::sqrt(s.real * s.real + s.imag * s.imag);
        case PlotMode::Phase: return std::atan2(s.real, s.imag);
    }
    return s.real;
}

PlotSession::PlotSession(const std::string& title, int line_thickness)
    : title_(title), line_thickness_(line_thickness) {}

void PlotSession::add_trace(const Trace& trace) {
    traces_.push_back(trace);
}

void PlotSession::set_axis_mode(XAxisMode mode) {
    xaxis_mode_ = mode;
}

void PlotSession::set_plot_mode(PlotMode mode) {
    mode_ = mode;
}

void PlotSession::set_zoom_range(ZoomRegion zr) {
    current_view_ = zr;
}

void PlotSession::set_fixed_y_range(std::optional<std::pair<double, double>> range) {
    fixed_y_range_ = range;
}

void PlotSession::set_fixed_x_range(std::optional<std::pair<double, double>> range) {
    fixed_x_range_ = range;
}

ZoomRegion PlotSession::autoscale_region(const std::vector<PlotSample>& samples) {
    double xmin, xmax;
    if (fixed_x_range_) {
        xmin = fixed_x_range_->first;
        xmax = fixed_x_range_->second;
    } else {
        if (mode_ == PlotMode::IQ){
            for (const auto& s : samples) {
                xmin = std::min(xmin, s.real);
                xmax = std::max(xmax, s.real);
            }
            if (xmax - xmin == 0) {
                xmin -= 1.0; xmax += 1.0;
            } else {
                xmin -= 0.05 * (xmax - xmin);
                xmax += 0.05 * (xmax - xmin);
            }
        } else {
            xmin = samples.front().time;
            xmax = samples.back().time;
        }
    }

    double ymin = compute_value(samples.front(), mode_);
    double ymax = ymin;

    if (fixed_y_range_) {
        ymin = fixed_y_range_->first;
        ymax = fixed_y_range_->second;
    } else {
        if (mode_ == PlotMode::IQ){
            for (const auto& s : samples) {
                ymin = std::min(ymin, s.imag);
                ymax = std::max(ymax, s.imag);
            }
        } else {
            for (const auto& s : samples) {
                double v = compute_value(s, mode_);
                ymin = std::min(ymin, v);
                ymax = std::max(ymax, v);
            }
        }
        if (ymax - ymin == 0) {
            ymin -= 1.0; ymax += 1.0;
        } else {
            ymin -= 0.05 * (ymax - ymin);
            ymax += 0.05 * (ymax - ymin);
        }
    }

    return {xmin, xmax, ymin, ymax};
}

void DecimatedTrace::update_if_needed(const std::vector<PlotSample>& samples,
                                      PlotMode mode,
                                      const ZoomRegion& view,
                                      int pixel_width) {
    if (!dirty &&
        cached_region.xmin == view.xmin &&
        cached_region.xmax == view.xmax &&
        cached_mode == mode &&
        cached_width == pixel_width) return;

    y_minmax_per_pixel.assign(pixel_width, {std::numeric_limits<double>::max(),
                                            std::numeric_limits<double>::lowest()});

    for (const auto& s : samples) {
        double xval = s.time;
        if (xval < view.xmin || xval > view.xmax) continue;

        int bin = static_cast<int>((xval - view.xmin) / (view.xmax - view.xmin) * pixel_width);
        if (bin < 0 || bin >= pixel_width) continue;

        double val = compute_value(s, mode);
        auto& [minv, maxv] = y_minmax_per_pixel[bin];
        minv = std::min(minv, val);
        maxv = std::max(maxv, val);
    }

    cached_mode = mode;
    cached_region = view;
    cached_width = pixel_width;
    dirty = false;
}

void PlotSession::save_pixmap_to_ppm(const char* filename) {
    XImage* image = XGetImage(dpy_, pixmap_, 0, 0, width_, height_, AllPlanes, ZPixmap);
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6\n" << width_ << " " << height_ << "\n255\n";
    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            unsigned long p = XGetPixel(image, x, y);
            ofs.put((p >> 16) & 0xFF);
            ofs.put((p >> 8) & 0xFF);
            ofs.put(p & 0xFF);
        }
    }
    XDestroyImage(image);
    std::cout << "Saved image to " << filename << "\\n";
}

void PlotSession::init_x11() {
    dpy_ = XOpenDisplay(nullptr);
    if (!dpy_) {
        std::cerr << "Unable to open X11 display\\n";
        exit(1);
    }

    screen_ = DefaultScreen(dpy_);
    win_ = XCreateSimpleWindow(dpy_, RootWindow(dpy_, screen_), 10, 10, width_, height_, 1,
                               BlackPixel(dpy_, screen_), WhitePixel(dpy_, screen_));
    XSelectInput(dpy_, win_, ExposureMask | ButtonPressMask | ButtonReleaseMask |
                          PointerMotionMask | StructureNotifyMask);
    XMapWindow(dpy_, win_);
    gc_ = XCreateGC(dpy_, win_, 0, nullptr);
    pixmap_ = XCreatePixmap(dpy_, win_, width_, height_, DefaultDepth(dpy_, screen_));
}

void PlotSession::draw_toolbar() {
    int tool_y = height_ - TOOLBAR_HEIGHT, xpos = 20;
    toolbar_buttons_ = {
        {"Save PNG", 0, 100},
        {"Cycle X-Axis", 0, 120},
        {"Magnitude", 0, 100},
        {"Real", 0, 100},
        {"Imag", 0, 100},
        {"Phase", 0, 100},
        {"Imag vs Real", 0, 120}
    };

    for (auto& b : toolbar_buttons_) {
        XDrawRectangle(dpy_, pixmap_, gc_, xpos, tool_y, 100, TOOLBAR_HEIGHT - 10);
        XDrawString(dpy_, pixmap_, gc_, xpos + 10, tool_y + 18, b.label.c_str(), b.label.size());
        b.x = xpos;
        b.width = 100;
        xpos += 120;
    }
}

void PlotSession::draw_legend() {
    legend_boxes_.clear();
    int plot_width = width_ * (1.0 - LEFT_MARGIN - RIGHT_MARGIN);
    int px = width_ * LEFT_MARGIN;
    int legend_start_x = px + 10;
    int legend_y = height_ - TOOLBAR_HEIGHT - VERTICAL_PADDING - LEGEND_HEIGHT;
    int box_height = 14;
    int box_width = 12;
    int spacing_x = 100;

    for (size_t t = 0; t < traces_.size(); ++t) {
        int lx = legend_start_x + t * spacing_x;

        if (!traces_[t].visible) {
            XSetForeground(dpy_, gc_, 0x404040);  // dark gray
        } else {
            XSetForeground(dpy_, gc_, trace_colors[t % trace_colors.size()]);
        }

        XFillRectangle(dpy_, pixmap_, gc_, lx, legend_y, box_width, box_height);

        XSetForeground(dpy_, gc_, COLOR_FG);
        std::string label = traces_[t].label;
        XDrawString(dpy_, pixmap_, gc_, lx + box_width + 6, legend_y + 10,
                    label.c_str(), label.size());

        legend_boxes_.push_back({lx, legend_y, spacing_x, box_height});
    }
}

void PlotSession::render_pixmap() {
    XSetForeground(dpy_, gc_, COLOR_BG);
    XFillRectangle(dpy_, pixmap_, gc_, 0, 0, width_, height_);

    int plot_width = width_ * (1.0 - LEFT_MARGIN - RIGHT_MARGIN);
    int plot_height = height_ * (1.0 - TOP_MARGIN - BOTTOM_MARGIN);
    int px = width_ * LEFT_MARGIN;
    int py = height_ * TOP_MARGIN;

    plot_height -= (TOOLBAR_HEIGHT + VERTICAL_PADDING + LEGEND_HEIGHT);

    if (mode_ == PlotMode::IQ) {
        // Draw axes (optional: add ticks/labels)
        XSetForeground(dpy_, gc_, COLOR_FG);
        XDrawLine(dpy_, pixmap_, gc_, px, py + plot_height / 2, px + plot_width, py + plot_height / 2); // x-axis
        XDrawLine(dpy_, pixmap_, gc_, px + plot_width / 2, py, px + plot_width / 2, py + plot_height); // y-axis

        // Draw traces as points in I/Q plane
        for (size_t t = 0; t < traces_.size(); ++t) {
            const auto& trace = traces_[t];
            XSetForeground(dpy_, gc_, trace_colors[t % trace_colors.size()]);
            for (const auto& s : trace.samples) {
                // Map real/imag to pixel coordinates
                int x = px + static_cast<int>(
                    (s.real - current_view_.xmin) / (current_view_.xmax - current_view_.xmin) * plot_width
                );
                int y = py + static_cast<int>(
                    (current_view_.ymax - s.imag) / (current_view_.ymax - current_view_.ymin) * plot_height
                );
                XFillArc(dpy_, pixmap_, gc_, x - line_thickness_, y - line_thickness_,
                         2 * line_thickness_, 2 * line_thickness_, 0, 360 * 64);
            }
        }
        // Optionally: draw legend, toolbar, etc.
        draw_legend();
        draw_toolbar();
        return;
    }

    if (decimated_cache_.size() != traces_.size()) {
        decimated_cache_.resize(traces_.size());
    }

    for (size_t i = 0; i < traces_.size(); ++i) {
        decimated_cache_[i].update_if_needed(traces_[i].samples, mode_, current_view_, plot_width);
    }

    XGCValues old_vals;
    XGetGCValues(dpy_, gc_, GCLineWidth | GCLineStyle | GCCapStyle | GCJoinStyle, &old_vals);
    XSetLineAttributes(dpy_, gc_, line_thickness_, LineSolid, CapButt, JoinMiter);

    for (size_t t = 0; t < traces_.size(); ++t) {
        if (!traces_[t].visible) continue;
        const auto& bins = decimated_cache_[t].y_minmax_per_pixel;
        if (bins.empty()) continue;

        XSetForeground(dpy_, gc_, trace_colors[t % trace_colors.size()]);
        int prev_x = -1, prev_y = -1;

        for (int xpix = 0; xpix < static_cast<int>(bins.size()); ++xpix) {
            const auto& [minv, maxv] = bins[xpix];
            if (minv > maxv) continue;

            int screen_x = px + xpix;
            int y1 = py + static_cast<int>((current_view_.ymax - minv) / (current_view_.ymax - current_view_.ymin) * plot_height);
            int y2 = py + static_cast<int>((current_view_.ymax - maxv) / (current_view_.ymax - current_view_.ymin) * plot_height);
            y1 = std::clamp(y1, py, py + plot_height - 1);
            y2 = std::clamp(y2, py, py + plot_height - 1);

            if (traces_[t].style == PlotRenderStyle::LINES) {
                XDrawLine(dpy_, pixmap_, gc_, screen_x, y1, screen_x, y2);
                if (prev_x != -1 && std::abs(screen_x - prev_x) >= 1) {
                    XDrawLine(dpy_, pixmap_, gc_, prev_x, prev_y, screen_x, y1);
                }
                prev_x = screen_x;
                prev_y = y2;
            } else {
                XFillArc(dpy_, pixmap_, gc_, screen_x - line_thickness_, y1 - line_thickness_,
                         2 * line_thickness_, 2 * line_thickness_, 0, 360 * 64);
                if (minv != maxv) {
                    XFillArc(dpy_, pixmap_, gc_, screen_x - line_thickness_, y2 - line_thickness_,
                             2 * line_thickness_, 2 * line_thickness_, 0, 360 * 64);
                }
            }
        }
    }

    XSetForeground(dpy_, gc_, COLOR_FG);
    XSetLineAttributes(dpy_, gc_, 1, LineSolid, CapButt, JoinMiter);
    XDrawRectangle(dpy_, pixmap_, gc_, px, py, plot_width, plot_height);
    draw_axes();

    int title_x = width_ / 2 - (title_.length() * 3);
    XDrawString(dpy_, pixmap_, gc_, title_x, py - 10, title_.c_str(), title_.size());

    draw_toolbar();
    draw_legend();
}

void PlotSession::draw_axes() {
    int plot_width = width_ * (1.0 - LEFT_MARGIN - RIGHT_MARGIN);
    int plot_height = height_ * (1.0 - TOP_MARGIN - BOTTOM_MARGIN);
    int px = width_ * LEFT_MARGIN;
    int py = height_ * TOP_MARGIN;

    plot_height -= (TOOLBAR_HEIGHT + VERTICAL_PADDING + LEGEND_HEIGHT);

    const int ticks = 5;

    for (int i = 0; i <= ticks; ++i) {
        double val = current_view_.xmin + i * (current_view_.xmax - current_view_.xmin) / ticks;
        int x = px + i * plot_width / ticks;

        XDrawLine(dpy_, pixmap_, gc_, x, py + plot_height, x, py + plot_height + 5);

        char label[32];
        if (xaxis_mode_ == XAxisMode::INDEX) {
            const auto& s0 = traces_.empty() ? PlotSample{0, 0, 0} : traces_[0].samples[0];
            double dt = traces_[0].samples[1].time - s0.time;
            int idx = static_cast<int>(std::round((val - s0.time) / dt));
            snprintf(label, sizeof(label), "%d", idx);
        } else {
            snprintf(label, sizeof(label), "%.2f", val);
        }

        XDrawString(dpy_, pixmap_, gc_, x - 10, py + plot_height + 20, label, strlen(label));
    }

    for (int i = 0; i <= ticks; ++i) {
        double v = current_view_.ymin + i * (current_view_.ymax - current_view_.ymin) / ticks;
        int y = py + static_cast<int>((current_view_.ymax - v) / (current_view_.ymax - current_view_.ymin) * plot_height);

        XDrawLine(dpy_, pixmap_, gc_, px - 5, y, px, y);

        char label[32];
        snprintf(label, sizeof(label), "%.2f", v);
        int text_width = XTextWidth(XQueryFont(dpy_, XGContextFromGC(gc_)), label, strlen(label));
        XDrawString(dpy_, pixmap_, gc_, px - text_width - 10, y + 4, label, strlen(label));
    }
}

void PlotSession::reset_zoom() {
    while (!zoom_stack_.empty()) zoom_stack_.pop();
    zoom_depth_ = 0;
    current_view_ = autoscale_region(traces_[0].samples);
}

void PlotSession::run() {
    if (traces_.empty()) return;

    init_x11();
    current_view_ = autoscale_region(traces_[0].samples);

    XEvent e;
    while (true) {
        if (XPending(dpy_)) {
            XNextEvent(dpy_, &e);

            if (e.type == ConfigureNotify) {
                if (e.xconfigure.width != width_ || e.xconfigure.height != height_) {
                    width_ = e.xconfigure.width;
                    height_ = e.xconfigure.height;
                    XFreePixmap(dpy_, pixmap_);
                    pixmap_ = XCreatePixmap(dpy_, win_, width_, height_, DefaultDepth(dpy_, screen_));
                    pixmap_dirty_ = true;
                }

            } else if (e.type == ButtonPress) {
                int x = e.xbutton.x;
                int y = e.xbutton.y;

                if (y >= height_ - TOOLBAR_HEIGHT - 10) {
                    for (auto& b : toolbar_buttons_) {
                        if (x >= b.x && x <= b.x + b.width) {
                            if (b.label == "Save PNG") {
                                save_pixmap_to_ppm("plot_out.ppm");
                            } else if (b.label == "Magnitude") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Magnitude;
                            } else if (b.label == "Imag") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Real;
                            } else if (b.label == "Real") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Imag;
                            } else if (b.label == "Phase") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Phase;
                            } else if (b.label == "Imag vs Real") {
                                if (mode_ != PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::IQ;
                            } else if (b.label == "Cycle X-Axis") {
                                xaxis_mode_ = (xaxis_mode_ == XAxisMode::TIME) ? XAxisMode::INDEX : XAxisMode::TIME;
                            }
                            pixmap_dirty_ = true;
                        }
                    }

                } else if (y >= height_ - TOOLBAR_HEIGHT - VERTICAL_PADDING - LEGEND_HEIGHT) {
                    for (size_t t = 0; t < legend_boxes_.size(); ++t) {
                        const auto& box = legend_boxes_[t];
                        if (y >= box.y && y <= box.y + box.height &&
                            x >= box.x && x <= box.x + box.width) {
                            if (e.xbutton.button == Button1)
                                traces_[t].visible = !traces_[t].visible;
                            else if (e.xbutton.button == Button3)
                                traces_[t].style = (traces_[t].style == PlotRenderStyle::LINES)
                                                    ? PlotRenderStyle::DOTS
                                                    : PlotRenderStyle::LINES;
                            pixmap_dirty_ = true;
                        }
                    }

                } else {
                    if (e.xbutton.button == Button1) {
                        dragging_ = true;
                        bx0_ = bx1_ = e.xbutton.x;
                        by0_ = by1_ = e.xbutton.y;
                    } else if (e.xbutton.button == Button3) {
                        if (!zoom_stack_.empty()) {
                            current_view_ = zoom_stack_.top();
                            zoom_stack_.pop();
                            zoom_depth_--;
                            pixmap_dirty_ = true;
                        }
                    }
                }

            } else if (e.type == ButtonRelease && dragging_) {
                dragging_ = false;
                bx1_ = e.xbutton.x;
                by1_ = e.xbutton.y;

                int px = width_ * LEFT_MARGIN;
                int py = height_ * TOP_MARGIN;
                int pw = width_ * (1.0 - LEFT_MARGIN - RIGHT_MARGIN);
                int ph = height_ * (1.0 - TOP_MARGIN - BOTTOM_MARGIN) - (TOOLBAR_HEIGHT + VERTICAL_PADDING + LEGEND_HEIGHT);

                if (zoom_depth_ < MAX_ZOOM_HISTORY) {
                    int x1 = std::min(bx0_, bx1_), x2 = std::max(bx0_, bx1_);
                    int y1 = std::min(by0_, by1_), y2 = std::max(by0_, by1_);

                    if (x2 - x1 > 10 && y2 - y1 > 10) {
                        zoom_stack_.push(current_view_);
                        zoom_depth_++;

                        double nxmin = current_view_.xmin + (x1 - px) / static_cast<double>(pw) * (current_view_.xmax - current_view_.xmin);
                        double nxmax = current_view_.xmin + (x2 - px) / static_cast<double>(pw) * (current_view_.xmax - current_view_.xmin);
                        double nymin = current_view_.ymax - (y2 - py) / static_cast<double>(ph) * (current_view_.ymax - current_view_.ymin);
                        double nymax = current_view_.ymax - (y1 - py) / static_cast<double>(ph) * (current_view_.ymax - current_view_.ymin);
                        current_view_ = {nxmin, nxmax, nymin, nymax};
                        pixmap_dirty_ = true;
                    }
                }

                bx0_ = by0_ = bx1_ = by1_ = -1;

            } else if (e.type == MotionNotify) {
                int mx = e.xmotion.x;
                int my = e.xmotion.y;

                if (dragging_) {
                    bx1_ = mx;
                    by1_ = my;
                }

                int plot_width = width_ * (1.0 - LEFT_MARGIN - RIGHT_MARGIN);
                int plot_height = height_ * (1.0 - TOP_MARGIN - BOTTOM_MARGIN);
                int px = width_ * LEFT_MARGIN;
                int py = height_ * TOP_MARGIN;

                plot_height -= (TOOLBAR_HEIGHT + VERTICAL_PADDING + LEGEND_HEIGHT);

                int rel_x = std::clamp(mx - px, 0, plot_width - 1);
                int rel_y = std::clamp(my - py, 0, plot_height - 1);

                double t = current_view_.xmin + (rel_x / (double)plot_width) * (current_view_.xmax - current_view_.xmin);
                double v = current_view_.ymax - (rel_y / (double)plot_height) * (current_view_.ymax - current_view_.ymin);

                std::ostringstream ss;
                ss.precision(2);
                if (xaxis_mode_ == XAxisMode::TIME) {
                    ss << "Time: " << t;
                } else {
                    double dt = traces_[0].samples[1].time - traces_[0].samples[0].time;
                    int idx = static_cast<int>(std::round((t - traces_[0].samples[0].time) / dt));
                    ss << "Index: " << idx;
                }
                ss << "  Value: " << v;
                readout_ = ss.str();
            }

        } else {
            if (pixmap_dirty_) {
                render_pixmap();
                pixmap_dirty_ = false;
            }

            XCopyArea(dpy_, pixmap_, win_, gc_, 0, 0, width_, height_, 0, 0);

            int label_y = height_ - TOOLBAR_HEIGHT + 10;
            std::string label_text = "Mode: " + mode_to_string(mode_) + "  " + readout_;
            XSetForeground(dpy_, gc_, COLOR_FG);
            XDrawString(dpy_, win_, gc_, width_ - READOUT_POS_FROM_RIGHT, label_y, label_text.c_str(), label_text.size());

            if (dragging_ && bx0_ >= 0 && bx1_ >= 0) {
                int zx = std::min(bx0_, bx1_);
                int zy = std::min(by0_, by1_);
                int zw = std::abs(bx1_ - bx0_);
                int zh = std::abs(by1_ - by0_);
                XSetForeground(dpy_, gc_, COLOR_FG);
                XDrawRectangle(dpy_, win_, gc_, zx, zy, zw, zh);
            }

            XFlush(dpy_);
            usleep(10000);
        }
    }

    XFreePixmap(dpy_, pixmap_);
    XFreeGC(dpy_, gc_);
    XDestroyWindow(dpy_, win_);
    XCloseDisplay(dpy_);
}

void PlotSession::render_once() {
    init_x11();
    current_view_ = autoscale_region(traces_[0].samples);
    render_pixmap();
    XCopyArea(dpy_, pixmap_, win_, gc_, 0, 0, width_, height_, 0, 0);
    XFlush(dpy_);
}

} // namespace xplot