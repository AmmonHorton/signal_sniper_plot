// Needs an X display (e.g. `Xvfb :99 &` then DISPLAY=:99). Skips when none is available.
//   DISPLAY=:99 bazel test //ssp/platform/x11:x11_backend_test --test_env=DISPLAY

// gtest and our headers first: X11 #defines names like None and Expose.
#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <thread>

#include "ssp/platform/backend.h"

#include <X11/Xlib.h>
#include <X11/keysym.h>

namespace ssp {
namespace {

using T = InputEvent::Type;

/// Second client connection that finds our window and sends it synthetic events.
class Poker {
public:
    Poker() : dpy_(XOpenDisplay(nullptr)) {}
    ~Poker() {
        if (dpy_) XCloseDisplay(dpy_);
    }
    Window find(const char* name) {
        Window root_ret, parent, *kids = nullptr;
        unsigned n = 0;
        XQueryTree(dpy_, DefaultRootWindow(dpy_), &root_ret, &parent, &kids, &n);
        Window found = 0;
        for (unsigned i = 0; i < n && !found; ++i) {
            char* wm = nullptr;
            if (XFetchName(dpy_, kids[i], &wm) && wm) {
                if (std::string(wm) == name) found = kids[i];
                XFree(wm);
            }
        }
        if (kids) XFree(kids);
        return found;
    }
    void key(Window w, KeySym sym, unsigned state = 0) {
        XEvent e{};
        e.xkey.type = KeyPress;
        e.xkey.window = w;
        e.xkey.display = dpy_;
        e.xkey.keycode = XKeysymToKeycode(dpy_, sym);
        e.xkey.state = state;
        e.xkey.x = 12;
        e.xkey.y = 34;
        XSendEvent(dpy_, w, True, KeyPressMask, &e);
        XFlush(dpy_);
    }
    void button(Window w, int type, unsigned b, int x, int y) {
        XEvent e{};
        e.xbutton.type = type;
        e.xbutton.window = w;
        e.xbutton.button = b;
        e.xbutton.x = x;
        e.xbutton.y = y;
        XSendEvent(dpy_, w, True, type == ButtonPress ? ButtonPressMask : ButtonReleaseMask, &e);
        XFlush(dpy_);
    }
    void close(Window w) {
        XEvent e{};
        e.xclient.type = ClientMessage;
        e.xclient.window = w;
        e.xclient.message_type = XInternAtom(dpy_, "WM_PROTOCOLS", False);
        e.xclient.format = 32;
        e.xclient.data.l[0] = static_cast<long>(XInternAtom(dpy_, "WM_DELETE_WINDOW", False));
        XSendEvent(dpy_, w, False, NoEventMask, &e);
        XFlush(dpy_);
    }
    Display* dpy_;
};

/// Collect events until one of `type` arrives (or ~2 s pass).
InputEvent wait_for(Backend& be, T type) {
    for (int i = 0; i < 100; ++i) {
        std::vector<InputEvent> ev;
        be.wait_events(ev, 20, -1);
        for (const auto& e : ev) {
            if (e.type == type) return e;
        }
    }
    ADD_FAILURE() << "event never arrived";
    return {};
}

TEST(X11Backend, TranslatesInputAndPresents) {
    if (!std::getenv("DISPLAY")) GTEST_SKIP() << "no DISPLAY";
    std::unique_ptr<Backend> be;
    try {
        be = open_x11(400, 300, "ssp-x11-test");
    } catch (const std::exception& e) {
        GTEST_SKIP() << e.what();
    }
    EXPECT_EQ(be->width(), 400);

    Framebuffer fb(be->width(), be->height(), 0x203040);
    fb.text(10, 10, "hello", 0xffffff);
    be->present(fb, {fb.bounds()});
    be->present(fb, {{5, 5, 20, 20}, {-50, -50, 10, 10}});  // partly/fully off-window rects are fine

    Poker poke;
    Window w = 0;
    for (int i = 0; i < 50 && !w; ++i) {
        w = poke.find("ssp-x11-test");
        if (!w) std::this_thread::sleep_for(std::chrono::milliseconds(20));
    }
    ASSERT_NE(w, 0u);

    poke.key(w, XK_Escape);
    InputEvent k = wait_for(*be, T::Key);
    EXPECT_EQ(k.key, key::kEscape);
    EXPECT_EQ(k.x, 12);

    poke.key(w, XK_s, ControlMask);
    k = wait_for(*be, T::Key);
    EXPECT_EQ(k.key, static_cast<uint32_t>('s'));
    EXPECT_EQ(k.mods, mods::kCtrl);

    poke.button(w, ButtonPress, 3, 100, 50);
    const InputEvent p = wait_for(*be, T::Press);
    EXPECT_EQ(p.button, 3);
    EXPECT_EQ(p.x, 100);
    EXPECT_EQ(p.y, 50);

    // Waiting with nothing pending honours the timeout.
    const auto t0 = std::chrono::steady_clock::now();
    std::vector<InputEvent> none;
    be->wait_events(none, 50, -1);
    EXPECT_GE(std::chrono::steady_clock::now() - t0, std::chrono::milliseconds(40));

    poke.close(w);
    wait_for(*be, T::Close);
}

}  // namespace
}  // namespace ssp
