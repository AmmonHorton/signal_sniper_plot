/// @file backend.h
/// @brief The only interface between the app and a window system (X11 today).
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "render/framebuffer.h"

namespace ssp {

/// @brief Key codes are X keysyms (identical to xkbcommon keysyms, so Wayland maps 1:1).
namespace key {
constexpr uint32_t kEscape = 0xff1b;
constexpr uint32_t kHome = 0xff50;
constexpr uint32_t kLeft = 0xff51;
constexpr uint32_t kUp = 0xff52;
constexpr uint32_t kRight = 0xff53;
constexpr uint32_t kDown = 0xff54;
constexpr uint32_t kSpace = 0x20;
}  // namespace key

namespace mods {
constexpr unsigned kShift = 1, kCtrl = 2, kAlt = 4;
}

struct InputEvent {
    enum class Type { Press, Release, Motion, Leave, Key, Resize, Exposed, Close };
    Type type = Type::Motion;
    int x = 0, y = 0;     ///< Pointer position (Press/Release/Motion/Key).
    int button = 0;       ///< 1 left, 2 middle, 3 right, 4/5 wheel.
    unsigned mods = 0;    ///< mods::k*
    uint32_t key = 0;     ///< Keysym (Key).
    int w = 0, h = 0;     ///< New window size (Resize).
    Rect rect;            ///< Exposed area.
};

class Backend {
public:
    virtual ~Backend() = default;

    virtual int width() const = 0;
    virtual int height() const = 0;

    /// @brief Append pending input to `out`, waiting up to `timeout_ms` (-1 = forever) for
    /// input or for `extra_fd` to become readable. Consecutive motions are coalesced.
    virtual void wait_events(std::vector<InputEvent>& out, int timeout_ms, int extra_fd) = 0;

    /// @brief Show the given areas of `fb` (window-sized) on screen.
    virtual void present(const Framebuffer& fb, const std::vector<Rect>& damage) = 0;
};

/// @throws std::runtime_error when no X display can be opened.
std::unique_ptr<Backend> open_x11(int width, int height, const std::string& title);

}  // namespace ssp
