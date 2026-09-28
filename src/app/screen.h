/// @file screen.h
/// @brief Everything one plot screen owns: data, settings, zoom stack, pointer state, and
/// the reduce result currently being shown. UI-thread only (the worker sees only a
/// ReduceResult and the content's pyramids).
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "app/settings.h"
#include "app/trace_content.h"
#include "render/frame.h"

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
    enum class Kind { XRange, YRange };
    Kind kind = Kind::XRange;
    std::string text;
    std::string error;
};

/// @brief Transient pointer/drag/popup state.
struct Interaction {
    int mx = -1, my = -1;  ///< Last pointer position.
    bool inside = false;   ///< Pointer is inside the window.
    bool dragging = false; ///< Left-drag zoom box in progress.
    int x0 = 0, y0 = 0;    ///< Drag start.
    std::optional<MenuState> menu;
    std::optional<PromptState> prompt;
    bool help = false;     ///< Keypress help overlay is showing.
};

struct Screen {
    Screen(std::unique_ptr<TraceContent> c, const PlotOptions& o);

    std::unique_ptr<TraceContent> content;
    PlotOptions opts;  ///< As given; home view is derived from these.
    Settings set;
    ViewStack views;
    Interaction ui;

    // What is on screen now.
    Layout layout;
    Range yshown{-1.0, 1.0};
    std::vector<Rect> legend_hits;          ///< Index = trace.
    std::shared_ptr<ReduceResult> result;   ///< Latest reduce (possibly still running).
    bool stopped = false;                   ///< User stopped the latest reduce early.
    std::string message;                    ///< One-line status, e.g. "saved foo.png".

    /// @brief Level-0 view for the current settings (fixed ranges only apply in their
    /// original coordinate system and mode).
    View home() const;
    void go_home() { views.reset(home()); }

    /// @brief y range to draw with `done` columns of the latest result finished.
    Range y_for(int done) const;

    /// @brief XView the next reduce needs for the current settings and layout.
    XView xview() const;

    int window_w() const { return layout.title.w; }
    int window_h() const { return layout.readout.bottom(); }

    /// @brief Pixel → data coordinates in the current view (x may be fractional index).
    double px_to_x(int px) const;
    double py_to_y(int py) const;
};

}  // namespace ssp
