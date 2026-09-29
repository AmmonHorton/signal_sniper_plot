/// @file frame.h
/// @brief Everything around the data: layout, axes, grid, title, legend. Shared by xplot and xraster.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "ssp/render/axis.h"
#include "ssp/render/framebuffer.h"
#include "ssp/types.h"

namespace ssp {

struct Theme {
    uint32_t bg = 0x000000;
    uint32_t fg = 0xFFFFFF;
    uint32_t grid = 0x505050;
    uint32_t dim = 0x606060;
    std::vector<uint32_t> palette = {0x2CFF05, 0xFF0000, 0x00C0FF, 0xFF00FF, 0xFFFF00,
                                     0xFFA500, 0x00FFFF, 0xFFFFFF, 0xFF69B4, 0xB22222};
    uint32_t trace_color(std::size_t i) const { return palette[i % palette.size()]; }
};

/// @brief Pixel regions of a plot window.
struct Layout {
    Rect plot;     ///< Data area.
    Rect title;    ///< Row above the data area.
    Rect xlabels;  ///< Row of x tick labels under the data area.
    Rect readout;  ///< Two text lines at the bottom of the window.
};

/// @brief Fixed margins (like SigPlot): the data area depends only on the window size, so
/// changing y labels never invalidates reduced columns. Tick labels are at most
/// kMaxTickLabelChars wide (see make_ticks).
/// @param colorbar leave room on the right for a raster colour bar.
Layout compute_layout(int window_w, int window_h, bool colorbar = false);

/// @brief Tick count that suits a data-area size in pixels.
int xdivisions(int plot_w);
int ydivisions(int plot_h);

/// @param ydown y increases downwards (rasters: frame 0 at the top).
void draw_axes(Framebuffer& fb, const Layout& l, double x0, double x1, const AxisTicks& xt,
               double y0, double y1, const AxisTicks& yt, bool grid, const Theme& th,
               bool ydown = false);

/// @brief Vertical colour bar right of the data area, low z at the bottom, with z ticks.
void draw_colorbar(Framebuffer& fb, const Layout& l, double z0, double z1,
                   const std::array<uint32_t, 256>& lut, const Theme& th);
void draw_title(Framebuffer& fb, const Layout& l, const std::string& title, const Theme& th);

struct LegendEntry {
    std::string name;
    uint32_t color = 0;
    Style style = Style::Lines;
    bool visible = true;
};
/// @brief Boxed legend in the top-right of the data area. Returns each entry's clickable rect.
std::vector<Rect> draw_legend(Framebuffer& fb, const Rect& plot,
                              const std::vector<LegendEntry>& entries, const Theme& th);

const char* cmode_name(CMode m);  ///< "Magnitude", "Real", ...

}  // namespace ssp
