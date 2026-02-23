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

constexpr unsigned long COLOR_BG = 0x000000; ///< Background colour (black).
constexpr unsigned long COLOR_FG = 0xFFFFFF; ///< Foreground colour (white).

constexpr double LEFT_MARGIN   = 0.08; ///< Fraction of window width reserved for y-axis labels.
constexpr double RIGHT_MARGIN  = 0.03; ///< Fraction of window width reserved on the right edge.
constexpr double TOP_MARGIN    = 0.05; ///< Fraction of window height reserved above the data area.

constexpr int TOOLBAR_HEIGHT    = 30; ///< Pixel height of the toolbar row.
constexpr int VERTICAL_PADDING  = 15; ///< Vertical pixel gap between adjacent UI rows.
constexpr int LEGEND_HEIGHT     = 20; ///< Pixel height of the legend row.
constexpr int AXIS_LABEL_HEIGHT = 20; ///< Pixels reserved below the data area for x-tick labels.
constexpr int MAX_ZOOM_HISTORY  =  5; ///< Maximum number of zoom levels that can be undone.
constexpr int STREAM_COLS_PER_FRAME = 50; ///< Pixel columns rendered per event-loop iteration during streaming.
constexpr int MIN_WIDTH         = 400; ///< Minimum window width enforced via WM size hints.
constexpr int MIN_HEIGHT        = 300; ///< Minimum window height enforced via WM size hints.

/// @name Toolbar geometry — shared by draw_toolbar() and compute_layout().
/// @{
constexpr int TOOLBAR_XSTART    =  20; ///< Left edge of the first toolbar button (pixels).
constexpr int TOOLBAR_BTN_W     = 100; ///< Pixel width of each toolbar button.
constexpr int TOOLBAR_BTN_STEP  = 120; ///< Pixel stride between the left edges of consecutive buttons.
constexpr int TOOLBAR_BTN_COUNT =   7; ///< Total number of toolbar buttons.
/// @}

static std::vector<unsigned long> trace_colors = {
    0x2CFF05, 0xFF0000, 0x00C0FF, 0xFF00FF, 0xFFFF00,
    0xFFA500, 0x00FFFF, 0xFFFFFF, 0xFF69B4, 0xB22222
};

std::string mode_to_string(PlotMode mode) {
    switch (mode) {
        case PlotMode::Magnitude: return "Magnitude";
        case PlotMode::Real:      return "Real";
        case PlotMode::Imag:      return "Imag";
        case PlotMode::Phase:     return "Phase";
        case PlotMode::IQ:        return "I/Q";
        default:                  return "Unknown";
    }
}

double compute_value(const PlotSample& s, PlotMode mode) {
    switch (mode) {
        case PlotMode::Real:      return s.real;
        case PlotMode::Imag:      return s.imag;
        case PlotMode::Magnitude: return std::sqrt(s.real * s.real + s.imag * s.imag);
        case PlotMode::Phase:     return std::atan2(s.real, s.imag);
    }
    return s.real;
}

// ---------------------------------------------------------------------------
// PlotSession construction / setters
// ---------------------------------------------------------------------------

PlotSession::PlotSession(const std::string& title, int line_thickness)
    : title_(title), line_thickness_(line_thickness) {}

void PlotSession::add_trace(const Trace& trace)  { traces_.push_back(trace); }
void PlotSession::set_axis_mode(XAxisMode mode)   { xaxis_mode_ = mode; }
void PlotSession::set_plot_mode(PlotMode mode)    { mode_ = mode; }
void PlotSession::set_zoom_range(ZoomRegion zr)   { current_view_ = zr; }

void PlotSession::set_fixed_y_range(std::optional<std::pair<double, double>> range) {
    fixed_y_range_ = range;
}

void PlotSession::set_fixed_x_range(std::optional<std::pair<double, double>> range) {
    fixed_x_range_ = range;
}

// ---------------------------------------------------------------------------
// Layout
// ---------------------------------------------------------------------------

PlotLayout PlotSession::compute_layout() const {
    PlotLayout l;

    // Bottom-up: toolbar → legend → axis labels → data area
    l.toolbar_y = height_ - TOOLBAR_HEIGHT;
    l.legend_y  = l.toolbar_y - VERTICAL_PADDING - LEGEND_HEIGHT;

    l.px = static_cast<int>(width_ * LEFT_MARGIN);
    l.py = static_cast<int>(height_ * TOP_MARGIN);

    int right_edge = width_ - static_cast<int>(width_ * RIGHT_MARGIN);
    l.pw = std::max(1, right_edge - l.px);

    // Data area bottom must leave room for axis tick labels before the legend
    l.ph = std::max(1, l.legend_y - VERTICAL_PADDING - AXIS_LABEL_HEIGHT - l.py);

    // Readout sits right of the last toolbar button
    l.readout_x = TOOLBAR_XSTART
                  + (TOOLBAR_BTN_COUNT - 1) * TOOLBAR_BTN_STEP
                  + TOOLBAR_BTN_W + 10;

    return l;
}

// ---------------------------------------------------------------------------
// Autoscale
// ---------------------------------------------------------------------------

ZoomRegion PlotSession::autoscale_region(const std::vector<PlotSample>& samples) {
    double xmin = 0.0, xmax = 0.0;
    if (fixed_x_range_) {
        xmin = fixed_x_range_->first;
        xmax = fixed_x_range_->second;
    } else {
        if (mode_ == PlotMode::IQ) {
            if (!samples.empty()) { xmin = samples.front().real; xmax = samples.front().real; }
            for (const auto& s : samples) {
                xmin = std::min(xmin, s.real);
                xmax = std::max(xmax, s.real);
            }
            if (xmax - xmin == 0) { xmin -= 1.0; xmax += 1.0; }
            else { xmin -= 0.05 * (xmax - xmin); xmax += 0.05 * (xmax - xmin); }
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
        if (mode_ == PlotMode::IQ) {
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
        if (ymax - ymin == 0) { ymin -= 1.0; ymax += 1.0; }
        else { ymin -= 0.05 * (ymax - ymin); ymax += 0.05 * (ymax - ymin); }
    }

    return {xmin, xmax, ymin, ymax};
}

// ---------------------------------------------------------------------------
// DecimatedTrace
// ---------------------------------------------------------------------------

void DecimatedTrace::update_if_needed(const std::vector<PlotSample>& samples,
                                      PlotMode mode,
                                      const ZoomRegion& view,
                                      int pixel_width) {
    if (!dirty &&
        cached_region.xmin == view.xmin &&
        cached_region.xmax == view.xmax &&
        cached_mode  == mode &&
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

    cached_mode   = mode;
    cached_region = view;
    cached_width  = pixel_width;
    dirty = false;
}

// ---------------------------------------------------------------------------
// PPM save
// ---------------------------------------------------------------------------

void PlotSession::save_pixmap_to_ppm(const char* filename) {
    XImage* image = XGetImage(dpy_, pixmap_, 0, 0, width_, height_, AllPlanes, ZPixmap);
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6\n" << width_ << " " << height_ << "\n255\n";
    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            unsigned long p = XGetPixel(image, x, y);
            ofs.put((p >> 16) & 0xFF);
            ofs.put((p >>  8) & 0xFF);
            ofs.put( p        & 0xFF);
        }
    }
    XDestroyImage(image);
    std::cout << "Saved image to " << filename << "\n";
}

// ---------------------------------------------------------------------------
// X11 init
// ---------------------------------------------------------------------------

void PlotSession::init_x11() {
    dpy_ = XOpenDisplay(nullptr);
    if (!dpy_) {
        std::cerr << "Unable to open X11 display\n";
        exit(1);
    }

    screen_ = DefaultScreen(dpy_);
    win_ = XCreateSimpleWindow(dpy_, RootWindow(dpy_, screen_), 10, 10, width_, height_, 1,
                               BlackPixel(dpy_, screen_), WhitePixel(dpy_, screen_));
    XSelectInput(dpy_, win_, ExposureMask | ButtonPressMask | ButtonReleaseMask |
                             PointerMotionMask | StructureNotifyMask);
    XMapWindow(dpy_, win_);
    XFlush(dpy_);   // make the blank window visible before any computation

    wm_delete_window_ = XInternAtom(dpy_, "WM_DELETE_WINDOW", False);
    XSetWMProtocols(dpy_, win_, &wm_delete_window_, 1);

    XSizeHints hints = {};
    hints.flags      = PMinSize;
    hints.min_width  = MIN_WIDTH;
    hints.min_height = MIN_HEIGHT;
    XSetWMNormalHints(dpy_, win_, &hints);

    gc_     = XCreateGC(dpy_, win_, 0, nullptr);
    pixmap_ = XCreatePixmap(dpy_, win_, width_, height_, DefaultDepth(dpy_, screen_));
}

// ---------------------------------------------------------------------------
// Draw helpers
// ---------------------------------------------------------------------------

void PlotSession::draw_toolbar(const PlotLayout& l) {
    toolbar_buttons_ = {
        {"Save PNG",     0, TOOLBAR_BTN_W},
        {"Cycle X-Axis", 0, TOOLBAR_BTN_W},
        {"Magnitude",    0, TOOLBAR_BTN_W},
        {"Real",         0, TOOLBAR_BTN_W},
        {"Imag",         0, TOOLBAR_BTN_W},
        {"Phase",        0, TOOLBAR_BTN_W},
        {"Imag vs Real", 0, TOOLBAR_BTN_W},
    };

    XSetLineAttributes(dpy_, gc_, 2, LineSolid, CapButt, JoinMiter);
    int xpos = TOOLBAR_XSTART;
    for (auto& b : toolbar_buttons_) {
        XDrawRectangle(dpy_, pixmap_, gc_, xpos, l.toolbar_y, TOOLBAR_BTN_W, TOOLBAR_HEIGHT - 10);
        XDrawString(dpy_, pixmap_, gc_, xpos + 10, l.toolbar_y + 13,
                    b.label.c_str(), b.label.size());
        b.x = xpos;
        xpos += TOOLBAR_BTN_STEP;
    }
    XSetLineAttributes(dpy_, gc_, 1, LineSolid, CapButt, JoinMiter);
}

void PlotSession::draw_legend(const PlotLayout& l) {
    legend_boxes_.clear();
    const int box_w   = 12;
    const int box_h   = 14;
    const int step_x  = 100;
    int lx0 = l.px + 10;

    for (size_t t = 0; t < traces_.size(); ++t) {
        int lx = lx0 + t * step_x;

        XSetForeground(dpy_, gc_, traces_[t].visible
                       ? trace_colors[t % trace_colors.size()]
                       : 0x404040u);
        XFillRectangle(dpy_, pixmap_, gc_, lx, l.legend_y, box_w, box_h);

        XSetForeground(dpy_, gc_, COLOR_FG);
        const std::string& label = traces_[t].label;
        XDrawString(dpy_, pixmap_, gc_, lx + box_w + 6, l.legend_y + 10,
                    label.c_str(), label.size());

        legend_boxes_.push_back({lx, l.legend_y, step_x, box_h});
    }
}

void PlotSession::draw_axes(const PlotLayout& l) {
    const int ticks = 5;

    // X ticks
    for (int i = 0; i <= ticks; ++i) {
        double val = current_view_.xmin + i * (current_view_.xmax - current_view_.xmin) / ticks;
        int x = l.px + i * l.pw / ticks;

        XDrawLine(dpy_, pixmap_, gc_, x, l.py + l.ph, x, l.py + l.ph + 5);

        char label[32];
        if (xaxis_mode_ == XAxisMode::INDEX) {
            const auto& s0 = traces_.empty() ? PlotSample{0, 0, 0} : traces_[0].samples[0];
            double dt = traces_[0].samples[1].time - s0.time;
            int idx = static_cast<int>(std::round((val - s0.time) / dt));
            snprintf(label, sizeof(label), "%d", idx);
        } else {
            snprintf(label, sizeof(label), "%.2f", val);
        }
        XDrawString(dpy_, pixmap_, gc_, x - 10, l.py + l.ph + AXIS_LABEL_HEIGHT,
                    label, strlen(label));
    }

    // Y ticks
    for (int i = 0; i <= ticks; ++i) {
        double v = current_view_.ymin + i * (current_view_.ymax - current_view_.ymin) / ticks;
        int y = l.py + static_cast<int>(
            (current_view_.ymax - v) / (current_view_.ymax - current_view_.ymin) * l.ph);

        XDrawLine(dpy_, pixmap_, gc_, l.px - 5, y, l.px, y);

        char label[32];
        snprintf(label, sizeof(label), "%.2f", v);
        XFontStruct* fi = XQueryFont(dpy_, XGContextFromGC(gc_));
        int tw = fi ? XTextWidth(fi, label, strlen(label)) : 0;
        if (fi) XFreeFontInfo(nullptr, fi, 1);
        XDrawString(dpy_, pixmap_, gc_, l.px - tw - 10, y + 4, label, strlen(label));
    }
}

// ---------------------------------------------------------------------------
// Render — seed one trace for streaming
// ---------------------------------------------------------------------------

void PlotSession::start_trace_stream(size_t t) {
    stream_active_trace_ = t;
    stream_xpix_         = 0;

    stream_trace_state_[t] = TraceRenderState{};
    auto& st = stream_trace_state_[t];

    const auto& l       = layout_;
    auto&       dc      = decimated_cache_[t];
    const auto& samples = traces_[t].samples;

    // Seed sample_cursor to first sample at or after xmin
    st.sample_cursor = static_cast<size_t>(
        std::lower_bound(samples.begin(), samples.end(), current_view_.xmin,
                         [](const PlotSample& s, double v) { return s.time < v; })
        - samples.begin());

    // For LINES: interpolate left-border prev position if the trace straddles xmin
    if (traces_[t].style == PlotRenderStyle::LINES && samples.size() >= 2) {
        const double y_range = current_view_.ymax - current_view_.ymin;
        if (y_range > 0.0 && st.sample_cursor > 0 && st.sample_cursor < samples.size()) {
            const PlotSample& a  = samples[st.sample_cursor - 1];
            const PlotSample& b  = samples[st.sample_cursor];
            const double      dt = b.time - a.time;
            const double y_left  = (std::abs(dt) < 1e-10)
                ? compute_value(a, mode_)
                : compute_value(a, mode_) +
                  std::clamp((current_view_.xmin - a.time) / dt, 0.0, 1.0) *
                  (compute_value(b, mode_) - compute_value(a, mode_));

            if (y_left >= current_view_.ymin && y_left <= current_view_.ymax) {
                st.prev_x      = l.px;
                st.prev_y      = std::clamp(
                    l.py + static_cast<int>((current_view_.ymax - y_left) / y_range * l.ph),
                    l.py, l.py + l.ph - 1);
                st.prev_data_y = y_left;
            }
        }
    }

    // If dirty: initialise bins to empty and record the view/mode/width they represent
    if (dc.dirty) {
        dc.y_minmax_per_pixel.assign(l.pw, {std::numeric_limits<double>::max(),
                                            std::numeric_limits<double>::lowest()});
        dc.cached_region = current_view_;
        dc.cached_mode   = mode_;
        dc.cached_width  = l.pw;
    }
}

// ---------------------------------------------------------------------------
// Render — init (clear + static elements + per-trace state seed)
// ---------------------------------------------------------------------------

void PlotSession::render_pixmap_init() {
    layout_ = compute_layout();
    const auto& l = layout_;

    stream_paused_    = false;
    streaming_active_ = false;

    XSetForeground(dpy_, gc_, COLOR_BG);
    XFillRectangle(dpy_, pixmap_, gc_, 0, 0, width_, height_);

    // ── IQ mode: render everything at once, no streaming ─────────────────────
    if (mode_ == PlotMode::IQ) {
        XSetForeground(dpy_, gc_, COLOR_FG);
        XDrawLine(dpy_, pixmap_, gc_,
                  l.px,            l.py + l.ph / 2,
                  l.px + l.pw,     l.py + l.ph / 2);
        XDrawLine(dpy_, pixmap_, gc_,
                  l.px + l.pw / 2, l.py,
                  l.px + l.pw / 2, l.py + l.ph);

        for (size_t t = 0; t < traces_.size(); ++t) {
            XSetForeground(dpy_, gc_, trace_colors[t % trace_colors.size()]);
            for (const auto& s : traces_[t].samples) {
                int x = l.px + static_cast<int>(
                    (s.real - current_view_.xmin) /
                    (current_view_.xmax - current_view_.xmin) * l.pw);
                int y = l.py + static_cast<int>(
                    (current_view_.ymax - s.imag) /
                    (current_view_.ymax - current_view_.ymin) * l.ph);
                XFillArc(dpy_, pixmap_, gc_,
                         x - line_thickness_, y - line_thickness_,
                         2 * line_thickness_, 2 * line_thickness_, 0, 360 * 64);
            }
        }
        draw_legend(l);
        draw_toolbar(l);
        XSetForeground(dpy_, gc_, COLOR_FG);
        int title_x = width_ / 2 - static_cast<int>(title_.length() * 3);
        XDrawString(dpy_, pixmap_, gc_, title_x, l.py - 10, title_.c_str(), title_.size());
        return;  // streaming_active_ stays false
    }

    // ── Cache invalidation: mark dirty for any trace whose view/mode/width changed ──
    if (decimated_cache_.size() != traces_.size())
        decimated_cache_.resize(traces_.size());
    for (size_t i = 0; i < traces_.size(); ++i) {
        auto& dc = decimated_cache_[i];
        if (dc.cached_region.xmin != current_view_.xmin ||
            dc.cached_region.xmax != current_view_.xmax ||
            dc.cached_mode         != mode_              ||
            dc.cached_width        != l.pw) {
            dc.dirty = true;
        }
    }

    // ── Draw static elements into the pixmap ──────────────────────────────────
    XSetForeground(dpy_, gc_, COLOR_FG);
    draw_axes(l);
    {
        int title_x = width_ / 2 - static_cast<int>(title_.length() * 3);
        XDrawString(dpy_, pixmap_, gc_, title_x, l.py - 10, title_.c_str(), title_.size());
    }
    draw_toolbar(l);
    draw_legend(l);

    // ── Initialise per-trace stream state and queue ───────────────────────────
    stream_trace_state_.assign(traces_.size(), TraceRenderState{});
    stream_pending_.clear();

    std::vector<size_t> visible;
    for (size_t t = 0; t < traces_.size(); ++t) {
        if (traces_[t].visible) visible.push_back(t);
    }

    if (!visible.empty()) {
        stream_pending_.assign(visible.begin() + 1, visible.end());
        start_trace_stream(visible[0]);
        streaming_active_ = true;
    }
}

// ---------------------------------------------------------------------------
// Render — stream next n pixel columns for all traces
// ---------------------------------------------------------------------------

void PlotSession::render_stream_columns(int n) {
    if (!streaming_active_) return;
    if (decimated_cache_.size() != traces_.size()) return;

    const auto& l = layout_;

    const size_t t = stream_active_trace_;
    if (t >= traces_.size()) { streaming_active_ = false; return; }

    auto&       dc      = decimated_cache_[t];
    auto&       st      = stream_trace_state_[t];
    const auto& samples = traces_[t].samples;

    // If the cache is already populated, render the whole trace in one frame
    const int effective_n = dc.dirty ? n : (l.pw + 1);
    const int col_end     = std::min(stream_xpix_ + effective_n, l.pw);

    const double y_range = current_view_.ymax - current_view_.ymin;
    const double xrange  = current_view_.xmax - current_view_.xmin;

    if (y_range <= 0.0 || xrange <= 0.0) { streaming_active_ = false; return; }

    XSetLineAttributes(dpy_, gc_, line_thickness_, LineSolid, CapButt, JoinMiter);
    XSetForeground(dpy_, gc_, trace_colors[t % trace_colors.size()]);

    for (int xpix = stream_xpix_; xpix < col_end; ++xpix) {
        double minv, maxv;

        if (dc.dirty) {
            // Advance sample_cursor: consume all samples whose bin is <= xpix.
            // A sample falls in bin floor((time-xmin)/xrange * pw), which is <= xpix
            // iff time < xmin + (xpix+1)*xrange/pw  (== col_xmax).
            const double col_xmax = current_view_.xmin + (xpix + 1) * xrange / l.pw;
            while (st.sample_cursor < samples.size() &&
                   samples[st.sample_cursor].time < col_xmax) {
                const double xval = samples[st.sample_cursor].time;
                if (xval >= current_view_.xmin && xval <= current_view_.xmax) {
                    const int bin = static_cast<int>(
                        (xval - current_view_.xmin) / xrange * l.pw);
                    if (bin >= 0 && bin < l.pw) {
                        const double val = compute_value(samples[st.sample_cursor], mode_);
                        auto& [bmin, bmax] = dc.y_minmax_per_pixel[bin];
                        bmin = std::min(bmin, val);
                        bmax = std::max(bmax, val);
                    }
                }
                ++st.sample_cursor;
            }
            minv = dc.y_minmax_per_pixel[xpix].first;
            maxv = dc.y_minmax_per_pixel[xpix].second;
        } else {
            minv = dc.y_minmax_per_pixel[xpix].first;
            maxv = dc.y_minmax_per_pixel[xpix].second;
        }

        if (minv > maxv) continue;  // empty bin

        const int sx = l.px + xpix;

        // ── Out-of-range bin ─────────────────────────────────────────────────
        if (maxv < current_view_.ymin || minv > current_view_.ymax) {
            const bool   above    = (minv > current_view_.ymax);
            const double oor_val  = above ? minv : maxv;
            const double boundary = above ? current_view_.ymax : current_view_.ymin;

            if (st.prev_x != -1 && traces_[t].style == PlotRenderStyle::LINES) {
                const double dv = oor_val - st.prev_data_y;
                if (std::abs(dv) > 1e-10) {
                    const double frac = std::clamp((boundary - st.prev_data_y) / dv, 0.0, 1.0);
                    const int x_exit  = std::clamp(
                        static_cast<int>(st.prev_x + frac * (sx - st.prev_x)),
                        l.px, l.px + l.pw - 1);
                    const int y_exit  = std::clamp(
                        l.py + static_cast<int>((current_view_.ymax - boundary) / y_range * l.ph),
                        l.py, l.py + l.ph - 1);
                    XDrawLine(dpy_, pixmap_, gc_, st.prev_x, st.prev_y, x_exit, y_exit);
                }
            }
            st.last_oor_x   = sx;
            st.last_oor_val = oor_val;
            st.prev_x       = -1;
            continue;
        }

        // ── In-range bin ─────────────────────────────────────────────────────
        const double draw_min = std::max(minv, current_view_.ymin);
        const double draw_max = std::min(maxv, current_view_.ymax);

        const int y1 = std::clamp(
            l.py + static_cast<int>((current_view_.ymax - draw_min) / y_range * l.ph),
            l.py, l.py + l.ph - 1);
        const int y2 = std::clamp(
            l.py + static_cast<int>((current_view_.ymax - draw_max) / y_range * l.ph),
            l.py, l.py + l.ph - 1);

        if (traces_[t].style == PlotRenderStyle::LINES) {
            XDrawLine(dpy_, pixmap_, gc_, sx, y1, sx, y2);

            if (st.last_oor_x != -1) {
                const bool   above    = (st.last_oor_val > current_view_.ymax);
                const double boundary = above ? current_view_.ymax : current_view_.ymin;
                const double dv       = draw_min - st.last_oor_val;
                if (std::abs(dv) > 1e-10) {
                    const double frac   = std::clamp((boundary - st.last_oor_val) / dv, 0.0, 1.0);
                    const int x_entry   = std::clamp(
                        static_cast<int>(st.last_oor_x + frac * (sx - st.last_oor_x)),
                        l.px, l.px + l.pw - 1);
                    const int y_entry   = std::clamp(
                        l.py + static_cast<int>((current_view_.ymax - boundary) / y_range * l.ph),
                        l.py, l.py + l.ph - 1);
                    XDrawLine(dpy_, pixmap_, gc_, x_entry, y_entry, sx, y1);
                }
                st.last_oor_x = -1;
            } else if (st.prev_x != -1 && std::abs(sx - st.prev_x) >= 1) {
                XDrawLine(dpy_, pixmap_, gc_, st.prev_x, st.prev_y, sx, y1);
            }

            st.prev_x      = sx;
            st.prev_y      = y2;
            st.prev_data_y = draw_max;
        } else {
            XFillArc(dpy_, pixmap_, gc_,
                     sx - line_thickness_, y1 - line_thickness_,
                     2 * line_thickness_, 2 * line_thickness_, 0, 360 * 64);
            if (draw_min != draw_max)
                XFillArc(dpy_, pixmap_, gc_,
                         sx - line_thickness_, y2 - line_thickness_,
                         2 * line_thickness_, 2 * line_thickness_, 0, 360 * 64);
            st.last_oor_x = -1;
        }
    }

    stream_xpix_ = col_end;

    // ── This trace is complete ────────────────────────────────────────────────
    if (stream_xpix_ >= l.pw) {
        // Draw right x-border exit segment for LINES
        if (st.prev_x != -1 && traces_[t].style == PlotRenderStyle::LINES && !samples.empty()) {
            auto it = std::lower_bound(samples.begin(), samples.end(), current_view_.xmax,
                                       [](const PlotSample& s, double v) { return s.time <= v; });
            if (it != samples.end() && it != samples.begin()) {
                const PlotSample& a = *(it - 1);
                const PlotSample& b = *it;
                const double dt = b.time - a.time;
                const double y_right = (std::abs(dt) < 1e-10)
                    ? compute_value(a, mode_)
                    : compute_value(a, mode_) +
                      std::clamp((current_view_.xmax - a.time) / dt, 0.0, 1.0) *
                      (compute_value(b, mode_) - compute_value(a, mode_));

                if (y_right >= current_view_.ymin && y_right <= current_view_.ymax) {
                    const int x_right    = l.px + l.pw - 1;
                    const int y_right_px = std::clamp(
                        l.py + static_cast<int>((current_view_.ymax - y_right) / y_range * l.ph),
                        l.py, l.py + l.ph - 1);
                    XDrawLine(dpy_, pixmap_, gc_, st.prev_x, st.prev_y, x_right, y_right_px);
                }
            }
        }

        dc.dirty = false;  // cache is now valid

        // Advance to the next pending trace, or finish
        if (!stream_pending_.empty()) {
            const size_t next = stream_pending_.front();
            stream_pending_.erase(stream_pending_.begin());
            start_trace_stream(next);
        } else {
            streaming_active_ = false;
        }
    }

    XSetLineAttributes(dpy_, gc_, 1, LineSolid, CapButt, JoinMiter);
}

// ---------------------------------------------------------------------------
// Zoom reset
// ---------------------------------------------------------------------------

void PlotSession::reset_zoom() {
    while (!zoom_stack_.empty()) zoom_stack_.pop();
    zoom_depth_ = 0;
    current_view_ = autoscale_region(traces_[0].samples);
}

// ---------------------------------------------------------------------------
// Event loop
// ---------------------------------------------------------------------------

void PlotSession::run() {
    if (traces_.empty()) return;

    init_x11();
    layout_ = compute_layout();
    current_view_ = autoscale_region(traces_[0].samples);

    XEvent e;
    while (true) {
        if (XPending(dpy_)) {
            XNextEvent(dpy_, &e);

            if (e.type == ClientMessage) {
                if (static_cast<Atom>(e.xclient.data.l[0]) == wm_delete_window_) break;

            } else if (e.type == ConfigureNotify) {
                if (e.xconfigure.width != width_ || e.xconfigure.height != height_) {
                    width_  = e.xconfigure.width;
                    height_ = e.xconfigure.height;
                    XFreePixmap(dpy_, pixmap_);
                    pixmap_ = XCreatePixmap(dpy_, win_, width_, height_,
                                            DefaultDepth(dpy_, screen_));
                    layout_ = compute_layout();
                    pixmap_dirty_ = true;
                }

            } else if (e.type == ButtonPress) {
                int x = e.xbutton.x;
                int y = e.xbutton.y;

                if (y >= layout_.toolbar_y) {
                    for (auto& b : toolbar_buttons_) {
                        if (x >= b.x && x <= b.x + b.width) {
                            if (b.label == "Save PNG") {
                                save_pixmap_to_ppm("plot_out.ppm");
                            } else if (b.label == "Magnitude") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Magnitude;
                            } else if (b.label == "Real") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Real;
                            } else if (b.label == "Imag") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Imag;
                            } else if (b.label == "Phase") {
                                if (mode_ == PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::Phase;
                            } else if (b.label == "Imag vs Real") {
                                if (mode_ != PlotMode::IQ) reset_zoom();
                                mode_ = PlotMode::IQ;
                            } else if (b.label == "Cycle X-Axis") {
                                xaxis_mode_ = (xaxis_mode_ == XAxisMode::TIME)
                                              ? XAxisMode::INDEX : XAxisMode::TIME;
                            }
                            pixmap_dirty_ = true;
                        }
                    }

                } else if (y >= layout_.legend_y) {
                    for (size_t t = 0; t < legend_boxes_.size(); ++t) {
                        const auto& box = legend_boxes_[t];
                        if (y >= box.y && y <= box.y + box.height &&
                            x >= box.x && x <= box.x + box.width) {
                            if (e.xbutton.button == Button1) {
                                traces_[t].visible = !traces_[t].visible;
                                if (!traces_[t].visible) {
                                    // Hide: pixels already in pixmap, need full redraw
                                    pixmap_dirty_ = true;
                                } else if (stream_trace_state_.size() == traces_.size() &&
                                           decimated_cache_.size() == traces_.size()) {
                                    // Show: additive — no pixmap clear needed
                                    const bool already_active =
                                        streaming_active_ && (stream_active_trace_ == t);
                                    const bool already_pending =
                                        std::find(stream_pending_.begin(),
                                                  stream_pending_.end(), t)
                                        != stream_pending_.end();
                                    if (!already_active && !already_pending) {
                                        if (!streaming_active_) {
                                            start_trace_stream(t);
                                            streaming_active_ = true;
                                        } else {
                                            stream_pending_.push_back(t);
                                        }
                                    }
                                    // Redraw legend so the color box lights back up
                                    XSetForeground(dpy_, gc_, COLOR_FG);
                                    draw_legend(layout_);
                                } else {
                                    // First render not yet completed — fall back
                                    pixmap_dirty_ = true;
                                }
                            } else if (e.xbutton.button == Button3) {
                                traces_[t].style = (traces_[t].style == PlotRenderStyle::LINES)
                                                   ? PlotRenderStyle::DOTS
                                                   : PlotRenderStyle::LINES;
                                pixmap_dirty_ = true;
                            }
                        }
                    }

                } else {
                    if (e.xbutton.button == Button1) {
                        dragging_ = true;
                        bx0_ = bx1_ = e.xbutton.x;
                        by0_ = by1_ = e.xbutton.y;
                    } else if (e.xbutton.button == Button2) {
                        stream_paused_ = !stream_paused_;
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

                if (zoom_depth_ < MAX_ZOOM_HISTORY) {
                    int x1 = std::min(bx0_, bx1_), x2 = std::max(bx0_, bx1_);
                    int y1 = std::min(by0_, by1_), y2 = std::max(by0_, by1_);

                    if (x2 - x1 > 3 && y2 - y1 > 3) {
                        zoom_stack_.push(current_view_);
                        zoom_depth_++;

                        double nxmin = current_view_.xmin + (x1 - layout_.px) /
                                       static_cast<double>(layout_.pw) *
                                       (current_view_.xmax - current_view_.xmin);
                        double nxmax = current_view_.xmin + (x2 - layout_.px) /
                                       static_cast<double>(layout_.pw) *
                                       (current_view_.xmax - current_view_.xmin);
                        double nymin = current_view_.ymax - (y2 - layout_.py) /
                                       static_cast<double>(layout_.ph) *
                                       (current_view_.ymax - current_view_.ymin);
                        double nymax = current_view_.ymax - (y1 - layout_.py) /
                                       static_cast<double>(layout_.ph) *
                                       (current_view_.ymax - current_view_.ymin);
                        current_view_ = {nxmin, nxmax, nymin, nymax};
                        pixmap_dirty_ = true;
                    }
                }
                bx0_ = by0_ = bx1_ = by1_ = -1;

            } else if (e.type == MotionNotify) {
                int mx = e.xmotion.x;
                int my = e.xmotion.y;

                mouse_x_ = mx;
                mouse_y_ = my;

                if (dragging_) { bx1_ = mx; by1_ = my; }

                int rel_x = std::clamp(mx - layout_.px, 0, layout_.pw - 1);
                int rel_y = std::clamp(my - layout_.py, 0, layout_.ph - 1);

                double t = current_view_.xmin +
                           (rel_x / static_cast<double>(layout_.pw)) *
                           (current_view_.xmax - current_view_.xmin);
                double v = current_view_.ymax -
                           (rel_y / static_cast<double>(layout_.ph)) *
                           (current_view_.ymax - current_view_.ymin);

                std::ostringstream ss;
                ss.precision(2);
                if (xaxis_mode_ == XAxisMode::TIME) {
                    ss << "Time: " << t;
                } else {
                    double dt = traces_[0].samples[1].time - traces_[0].samples[0].time;
                    int idx = static_cast<int>(
                        std::round((t - traces_[0].samples[0].time) / dt));
                    ss << "Index: " << idx;
                }
                ss << "  Value: " << v;
                readout_ = ss.str();
            }

        } else {
            if (pixmap_dirty_) {
                render_pixmap_init();
                pixmap_dirty_ = false;
            }

            if (streaming_active_ && !stream_paused_)
                render_stream_columns(STREAM_COLS_PER_FRAME);

            XCopyArea(dpy_, pixmap_, win_, gc_, 0, 0, width_, height_, 0, 0);

            // Plot border drawn on the window so it always sits on top of streaming data
            XSetForeground(dpy_, gc_, COLOR_FG);
            XSetLineAttributes(dpy_, gc_, 2, LineSolid, CapButt, JoinMiter);
            XDrawRectangle(dpy_, win_, gc_, layout_.px, layout_.py, layout_.pw, layout_.ph);
            XSetLineAttributes(dpy_, gc_, 1, LineSolid, CapButt, JoinMiter);

            std::string label_text = "Mode: " + mode_to_string(mode_) + "  " + readout_;
            XSetForeground(dpy_, gc_, COLOR_FG);
            XDrawString(dpy_, win_, gc_, layout_.readout_x, layout_.toolbar_y + 18,
                        label_text.c_str(), label_text.size());

            // Dotted crosshairs — only inside the plot area
            if (mouse_x_ >= layout_.px && mouse_x_ < layout_.px + layout_.pw &&
                mouse_y_ >= layout_.py && mouse_y_ < layout_.py + layout_.ph) {
                static const char dash_pattern[] = {4, 4};
                XSetLineAttributes(dpy_, gc_, 1, LineOnOffDash, CapButt, JoinMiter);
                XSetDashes(dpy_, gc_, 0, dash_pattern, 2);
                XSetForeground(dpy_, gc_, COLOR_FG);
                XDrawLine(dpy_, win_, gc_,
                          layout_.px,              mouse_y_,
                          layout_.px + layout_.pw, mouse_y_);
                XDrawLine(dpy_, win_, gc_,
                          mouse_x_, layout_.py,
                          mouse_x_, layout_.py + layout_.ph);
                XSetLineAttributes(dpy_, gc_, 1, LineSolid, CapButt, JoinMiter);
            }

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
    render_pixmap_init();
    while (streaming_active_)
        render_stream_columns(layout_.pw + 1);  // one call completes one trace
    XCopyArea(dpy_, pixmap_, win_, gc_, 0, 0, width_, height_, 0, 0);
    XSetForeground(dpy_, gc_, COLOR_FG);
    XSetLineAttributes(dpy_, gc_, 2, LineSolid, CapButt, JoinMiter);
    XDrawRectangle(dpy_, win_, gc_, layout_.px, layout_.py, layout_.pw, layout_.ph);
    XSetLineAttributes(dpy_, gc_, 1, LineSolid, CapButt, JoinMiter);
    XFlush(dpy_);
}

} // namespace xplot
