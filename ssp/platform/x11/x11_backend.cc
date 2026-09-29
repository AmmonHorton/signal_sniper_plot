// X11 (Xlib) backend: one window, client-side framebuffer uploaded with XPutImage.
//
// Only damaged rectangles are uploaded, so mouse-driven overlays stay cheap even over
// `ssh -X`. MIT-SHM could replace XPutImage later for faster full-frame uploads locally.

#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <poll.h>

#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include "ssp/platform/backend.h"

namespace ssp {
namespace {

constexpr int kMinW = 320, kMinH = 200;

/// Bit position of the lowest set bit of an 8-bit-per-channel mask.
int shift_of(unsigned long mask) {
    int s = 0;
    while (mask && !(mask & 1)) {
        mask >>= 1;
        ++s;
    }
    return s;
}

class X11Backend final : public Backend {
public:
    X11Backend(int w, int h, const std::string& title) : w_(std::max(w, kMinW)), h_(std::max(h, kMinH)) {
        dpy_ = XOpenDisplay(nullptr);
        if (!dpy_) throw std::runtime_error("cannot open X display (is DISPLAY set?)");
        const int scr = DefaultScreen(dpy_);
        visual_ = DefaultVisual(dpy_, scr);
        depth_ = DefaultDepth(dpy_, scr);
        if (visual_->c_class != TrueColor || depth_ < 24) {
            XCloseDisplay(dpy_);
            throw std::runtime_error("X display must be TrueColor with depth >= 24");
        }
        rs_ = shift_of(visual_->red_mask);
        gs_ = shift_of(visual_->green_mask);
        bs_ = shift_of(visual_->blue_mask);

        XSetWindowAttributes a{};
        a.background_pixmap = None;          // no server-side clears: no flashing on resize
        a.bit_gravity = NorthWestGravity;
        a.event_mask = ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask |
                       PointerMotionMask | LeaveWindowMask | StructureNotifyMask;
        win_ = XCreateWindow(dpy_, RootWindow(dpy_, scr), 0, 0, w_, h_, 0, depth_, InputOutput,
                             visual_, CWBackPixmap | CWBitGravity | CWEventMask, &a);
        gc_ = XCreateGC(dpy_, win_, 0, nullptr);

        XStoreName(dpy_, win_, title.c_str());
        XChangeProperty(dpy_, win_, XInternAtom(dpy_, "_NET_WM_NAME", False),
                        XInternAtom(dpy_, "UTF8_STRING", False), 8, PropModeReplace,
                        reinterpret_cast<const unsigned char*>(title.data()),
                        static_cast<int>(title.size()));
        XClassHint cls{const_cast<char*>("ssp"), const_cast<char*>("SignalSniperPlot")};
        XSetClassHint(dpy_, win_, &cls);
        XSizeHints hints{};
        hints.flags = PMinSize;
        hints.min_width = kMinW;
        hints.min_height = kMinH;
        XSetWMNormalHints(dpy_, win_, &hints);
        wm_delete_ = XInternAtom(dpy_, "WM_DELETE_WINDOW", False);
        XSetWMProtocols(dpy_, win_, &wm_delete_, 1);
        cursor_ = XCreateFontCursor(dpy_, XC_crosshair);
        XDefineCursor(dpy_, win_, cursor_);

        XMapWindow(dpy_, win_);
        XFlush(dpy_);
    }

    ~X11Backend() override {
        free_image();
        XFreeCursor(dpy_, cursor_);
        XFreeGC(dpy_, gc_);
        XDestroyWindow(dpy_, win_);
        XCloseDisplay(dpy_);
    }

    int width() const override { return w_; }
    int height() const override { return h_; }

    void wait_events(std::vector<InputEvent>& out, int timeout_ms, int extra_fd) override {
        if (!XPending(dpy_)) {
            pollfd fds[2] = {{ConnectionNumber(dpy_), POLLIN, 0}, {extra_fd, POLLIN, 0}};
            ::poll(fds, extra_fd >= 0 ? 2 : 1, timeout_ms);
        }
        while (XPending(dpy_)) {
            XEvent e;
            XNextEvent(dpy_, &e);
            translate(e, out);
        }
    }

    void present(const Framebuffer& fb, const std::vector<Rect>& damage) override {
        ensure_image(fb.width(), fb.height());
        for (Rect r : damage) {
            r = r.intersect(fb.bounds());
            if (r.empty()) continue;
            copy_to_image(fb, r);
            XPutImage(dpy_, win_, gc_, img_, r.x, r.y, r.x, r.y, r.w, r.h);
        }
        XFlush(dpy_);
    }

private:
    void translate(const XEvent& e, std::vector<InputEvent>& out) {
        using T = InputEvent::Type;
        InputEvent ev;
        switch (e.type) {
            case ButtonPress:
            case ButtonRelease:
                ev.type = e.type == ButtonPress ? T::Press : T::Release;
                ev.x = e.xbutton.x;
                ev.y = e.xbutton.y;
                ev.button = static_cast<int>(e.xbutton.button);
                ev.mods = mods_of(e.xbutton.state);
                break;
            case MotionNotify:
                ev.type = T::Motion;
                ev.x = e.xmotion.x;
                ev.y = e.xmotion.y;
                ev.mods = mods_of(e.xmotion.state);
                if (!out.empty() && out.back().type == T::Motion) {
                    out.back() = ev;  // coalesce
                    return;
                }
                break;
            case LeaveNotify:
                ev.type = T::Leave;
                break;
            case KeyPress: {
                XKeyEvent k = e.xkey;
                KeySym sym = 0;
                char buf[8];
                XLookupString(&k, buf, sizeof buf, &sym, nullptr);
                ev.type = T::Key;
                ev.key = static_cast<uint32_t>(sym);
                ev.x = k.x;
                ev.y = k.y;
                ev.mods = mods_of(k.state);
                break;
            }
            case ConfigureNotify:
                if (e.xconfigure.width == w_ && e.xconfigure.height == h_) return;
                w_ = e.xconfigure.width;
                h_ = e.xconfigure.height;
                ev.type = T::Resize;
                ev.w = w_;
                ev.h = h_;
                if (!out.empty() && out.back().type == T::Resize) {
                    out.back() = ev;
                    return;
                }
                break;
            case Expose:
                ev.type = T::Exposed;
                ev.rect = {e.xexpose.x, e.xexpose.y, e.xexpose.width, e.xexpose.height};
                break;
            case ClientMessage:
                if (static_cast<Atom>(e.xclient.data.l[0]) != wm_delete_) return;
                ev.type = T::Close;
                break;
            default:
                return;
        }
        out.push_back(ev);
    }

    static unsigned mods_of(unsigned state) {
        return ((state & ShiftMask) ? mods::kShift : 0) | ((state & ControlMask) ? mods::kCtrl : 0) |
               ((state & Mod1Mask) ? mods::kAlt : 0);
    }

    void ensure_image(int w, int h) {
        if (img_ && img_->width == w && img_->height == h) return;
        free_image();
        char* data = static_cast<char*>(std::malloc(static_cast<std::size_t>(w) * h * 4));
        img_ = XCreateImage(dpy_, visual_, depth_, ZPixmap, 0, data, w, h, 32, 0);
        if (!img_ || img_->bits_per_pixel != 32) {
            throw std::runtime_error("X11: only 32 bits-per-pixel images are supported");
        }
        // Straight copy when the server's pixel layout matches ours (0x00RRGGBB, host order).
        const bool host_lsb = [] { const uint16_t one = 1; return *reinterpret_cast<const uint8_t*>(&one) == 1; }();
        direct_ = rs_ == 16 && gs_ == 8 && bs_ == 0 && (img_->byte_order == LSBFirst) == host_lsb;
    }

    void free_image() {
        if (img_) XDestroyImage(img_);  // also frees the pixel data
        img_ = nullptr;
    }

    void copy_to_image(const Framebuffer& fb, const Rect& r) {
        for (int y = r.y; y < r.bottom(); ++y) {
            const uint32_t* src = fb.row(y) + r.x;
            char* dst = img_->data + static_cast<std::size_t>(y) * img_->bytes_per_line + r.x * 4;
            if (direct_) {
                std::memcpy(dst, src, static_cast<std::size_t>(r.w) * 4);
                continue;
            }
            for (int x = 0; x < r.w; ++x) {
                const uint32_t p = src[x];
                const unsigned long pix = ((p >> 16 & 0xff) << rs_) | ((p >> 8 & 0xff) << gs_) |
                                          ((p & 0xff) << bs_);
                XPutPixel(img_, r.x + x, y, pix);
            }
        }
    }

    Display* dpy_ = nullptr;
    Visual* visual_ = nullptr;
    int depth_ = 24;
    Window win_ = 0;
    GC gc_ = nullptr;
    Atom wm_delete_ = 0;
    Cursor cursor_ = 0;
    XImage* img_ = nullptr;
    bool direct_ = false;
    int rs_ = 16, gs_ = 8, bs_ = 0;
    int w_, h_;
};

}  // namespace

std::unique_ptr<Backend> open_x11(int width, int height, const std::string& title) {
    return std::make_unique<X11Backend>(width, height, title);
}

}  // namespace ssp
