/// @file controller.h
/// @brief Input → state changes. Pure with respect to the window system, so it is unit-tested
/// with synthetic events.
#pragma once

#include "app/screen.h"
#include "platform/backend.h"

namespace ssp {

/// @brief What the session has to do after an event. Bit flags.
enum Inval : unsigned {
    kNone = 0,
    kOverlay = 1 << 0,  ///< Crosshair, zoom box, readout or status changed.
    kRepaint = 1 << 1,  ///< Recompose from the existing reduce result.
    kReduce = 1 << 2,   ///< The view or plotted quantity changed: reduce again.
    kCancel = 1 << 3,   ///< Stop the running reduce, keep what is drawn.
    kResume = 1 << 4,   ///< Re-run the current view's reduce.
    kSave = 1 << 5,     ///< Save the plot as PNG.
    kQuit = 1 << 6,
};

/// @param rendering true while a reduce is still running.
unsigned handle_event(Screen& s, const InputEvent& e, bool rendering);

/// @brief Operations shared by keys and the menu. Each returns the Inval flags it needs.
namespace action {
unsigned set_mode(Screen& s, CMode m);
unsigned set_phunits(Screen& s, PhaseUnits u);
unsigned toggle_trace(Screen& s, std::size_t t);
unsigned toggle_index(Screen& s);
unsigned unzoom_all(Screen& s);
unsigned autoscale_y(Screen& s);
/// @brief Push `v` as a new zoom level (so right-click undoes it).
unsigned zoom_to(Screen& s, const View& v);
unsigned open_prompt(Screen& s, PromptState::Kind kind);
}  // namespace action

}  // namespace ssp
