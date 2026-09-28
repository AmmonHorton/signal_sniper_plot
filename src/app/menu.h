/// @file menu.h
/// @brief Middle-click popup menu: model (rebuilt from the screen's state each time), geometry,
/// input handling and drawing.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "app/screen.h"
#include "platform/backend.h"
#include "render/frame.h"

namespace ssp {

/// @brief One entry: a separator (empty label), a submenu (`sub` non-empty) or an action.
struct MenuItem {
    std::string label;
    std::string key;       ///< Shortcut shown on the right, e.g. "Home".
    bool checked = false;  ///< Shown with a '*' (current mode, enabled toggle).
    std::function<unsigned(Screen&)> action;
    std::vector<MenuItem> sub;

    bool separator() const { return label.empty(); }
};

std::vector<MenuItem> build_menu(const Screen& s);

/// @brief Open the menu with its corner at (x, y).
unsigned open_menu(Screen& s, int x, int y);

/// @brief Input while the menu is open; the menu consumes everything it is given.
unsigned handle_menu_event(Screen& s, const InputEvent& e);

/// @brief Draw the open menu; returns the rects it covered.
std::vector<Rect> draw_menu(Framebuffer& fb, const Screen& s, const Theme& th);

}  // namespace ssp
