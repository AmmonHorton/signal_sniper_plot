// xplot_main.cpp
// Minimal X-Midas xplot C++17 skeleton using Xlib

#include <X11/Xlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

// Minimal data layer (inspired by XPLOTLAYER)
struct XPlotLayer {
    std::vector<double> xbuf;
    std::vector<double> ybuf;
    bool display = true;
    int color = 0; // use XColor or RGB triplet later
};

// Global state (inspired by GX)
struct GXState {
    Display* display = nullptr;
    Window window;
    GC gc;
    int screen;
    std::vector<XPlotLayer> layers;
};

bool init_xwindow(GXState& gx, int width, int height) {
    gx.display = XOpenDisplay(NULL);
    if (!gx.display) {
        std::cerr << "Cannot open display" << std::endl;
        return false;
    }
    gx.screen = DefaultScreen(gx.display);
    gx.window = XCreateSimpleWindow(
        gx.display,
        RootWindow(gx.display, gx.screen),
        10, 10, width, height,
        1,
        BlackPixel(gx.display, gx.screen),
        WhitePixel(gx.display, gx.screen)
    );

    XSelectInput(gx.display, gx.window, ExposureMask | KeyPressMask);
    XMapWindow(gx.display, gx.window);

    gx.gc = XCreateGC(gx.display, gx.window, 0, NULL);
    return true;
}

void draw_plot(const GXState& gx, int width, int height) {
    for (const auto& layer : gx.layers) {
        if (!layer.display) continue;

        for (size_t i = 1; i < layer.xbuf.size(); ++i) {
            int x1 = static_cast<int>(layer.xbuf[i - 1] * width);
            int y1 = height - static_cast<int>(layer.ybuf[i - 1] * height);
            int x2 = static_cast<int>(layer.xbuf[i] * width);
            int y2 = height - static_cast<int>(layer.ybuf[i] * height);
            XDrawLine(gx.display, gx.window, gx.gc, x1, y1, x2, y2);
        }
    }
}

int main() {
    GXState gx;
    const int width = 800;
    const int height = 600;

    if (!init_xwindow(gx, width, height)) return EXIT_FAILURE;

    // Example data
    XPlotLayer layer;
    for (int i = 0; i < 100; ++i) {
        double x = i / 100.0;
        layer.xbuf.push_back(x);
        layer.ybuf.push_back(0.5 + 0.4 * std::sin(10 * x));
    }
    gx.layers.push_back(layer);

    // Main event loop
    XEvent e;
    while (true) {
        XNextEvent(gx.display, &e);
        if (e.type == Expose) {
            draw_plot(gx, width, height);
        } else if (e.type == KeyPress) {
            break;
        }
    }

    XDestroyWindow(gx.display, gx.window);
    XCloseDisplay(gx.display);
    return EXIT_SUCCESS;
}
