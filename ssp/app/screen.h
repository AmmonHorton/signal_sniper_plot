/// @file screen.h
/// @brief Everything one plot screen owns: data, settings, zoom stack, pointer state, and
/// the reduce result currently being shown. UI-thread only (the worker sees only a
/// ReduceResult and the content's pyramids).
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "ssp/app/content.h"
#include "ssp/app/settings.h"
#include "ssp/app/trace_content.h"
#include "ssp/render/frame.h"

namespace ssp {

/// @brief An open popup menu: its top-left corner and the highlighted item of each open
/// level (sel.size() levels are open; every level but the last highlights a submenu).
/// -1 means nothing highlighted.
struct MenuState {
    int x = 0, y = 0;
    std::vector<int> sel{-1};
};

/// @brief A one-line text prompt shown in the readout area.
struct PromptState {
    enum class Kind { XRange, YRange, ZRange };
    Kind kind = Kind::XRange;
    std::string text;
    std::string error;
};

class RasterContent;

/// @brief What a cut screen shows: one row (x-cut) or column (y-cut) of a raster, limited to
/// the raster's zoom box when the cut was made.
struct CutInfo {
    const RasterContent* source = nullptr;  ///< The raster screen below this one owns it.
    bool xcut = true;                       ///< Row (x-cut) or column (y-cut).
    std::size_t index = 0;                  ///< Which row / column.
    SampleRange steps;                      ///< Rows (x-cut) / columns (y-cut) in the box.
    SampleRange span;                       ///< Columns (x-cut) / rows (y-cut) in the box.
    bool index_axes = false;                ///< Raster was in index mode: label cells by number.
    std::string raster_title;
};

/// @brief Transient pointer/drag/popup state.
struct Interaction {
    int mx = -1, my = -1;  ///< Last pointer position.
    bool inside = false;   ///< Pointer is inside the window.
    bool dragging = false; ///< Left-drag zoom box in progress.
    int x0 = 0, y0 = 0;    ///< Drag start.
    std::optional<std::pair<double, double>> marker;  ///< Set by a left click (data coords).
    bool panning = false;  ///< Shift+left-drag pan in progress (from x0, y0).
    Range pan_x, pan_y;    ///< View when the pan started.
    std::optional<MenuState> menu;
    std::optional<PromptState> prompt;
    bool help = false;     ///< Keypress help overlay is showing.
};

struct Screen {
    explicit Screen(std::unique_ptr<Content> c);

    std::unique_ptr<Content> content;
    Settings set;
    ViewStack views;
    Interaction ui;

    // What is on screen now.
    std::vector<SampleRange> ir;  ///< IR: per-trace samples (the time window when IR was chosen).

    Layout layout;
    Range xshown{-1.0, 1.0};
    Range yshown{-1.0, 1.0};
    Range zshown{0.0, 1.0};  ///< Raster colour range on screen.
    std::vector<Rect> legend_hits;          ///< Index = trace.
    std::shared_ptr<ReduceResult> result;   ///< Latest reduce (possibly still running).
    bool stopped = false;                   ///< User stopped the latest reduce early.
    std::string message;                    ///< One-line status, e.g. "saved foo.png".

    std::optional<CutInfo> cut;  ///< Set on x/y-cut screens.
    /// Controller → session requests: open a cut of this row/column, or step the cut.
    std::optional<std::pair<bool, std::size_t>> cut_request;  ///< {xcut, index}
    int cut_step = 0;

    /// @brief Level-0 view for the current settings (fixed ranges only apply in their
    /// original coordinate system and mode).
    View home() const { return content->home(set); }
    /// @brief y increases downwards (rasters draw frame 0 at the top).
    bool ydown() const { return content->is_raster(); }
    Layout layout_for(int w, int h) const { return compute_layout(w, h, content->is_raster()); }
    /// @brief The xplot content, or null for a raster.
    TraceContent* traces() const { return dynamic_cast<TraceContent*>(content.get()); }
    /// @brief The raster content, or null for xplot.
    const RasterContent* raster() const;
    void go_home() { views.reset(home()); }

    /// @brief x/y range to draw with `done` columns of the latest result finished.
    Range x_for(int done) const;
    Range y_for(int done) const;
    Range z_for(int done) const;

    /// @brief XView the next reduce needs for the current settings and layout.
    XView xview() const;
    /// @brief Everything the next reduce for the top zoom level depends on.
    ReduceRequest request() const;
    /// @brief x range a level shows: its own, or for an autoscaled (IR home) level the range
    /// its finished result chose.
    Range level_x(const View& v) const;

    int window_w() const { return layout.title.w; }
    int window_h() const { return layout.readout.bottom(); }

    /// @brief Pixel → data coordinates in the current view (x may be fractional index).
    double px_to_x(int px) const;
    double py_to_y(int py) const;
};

}  // namespace ssp
