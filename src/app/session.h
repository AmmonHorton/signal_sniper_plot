/// @file session.h
/// @brief The interactive event loop, independent of the window system (tests drive it with a
/// fake Backend).
#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "platform/backend.h"
#include "ssp/plot.h"

namespace ssp {

class Content;

/// @brief Run an interactive session for `content` on `backend` until the window closes, `q`
/// is pressed, or `interrupt` (polled about every 100 ms when given) returns true.
void run_session(std::unique_ptr<Content> content, Backend& backend,
                 const std::function<bool()>& interrupt = {});

/// @brief Run an interactive xplot on `backend` until the window closes, `q` is pressed, or
/// `interrupt` (polled about every 100 ms when given) returns true.
void run_plot(std::vector<Signal> signals, const PlotOptions& options, Backend& backend,
              const std::function<bool()>& interrupt = {});

}  // namespace ssp
