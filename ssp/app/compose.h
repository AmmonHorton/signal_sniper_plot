/// @file compose.h
/// @brief Draw a whole plot (no overlays) into a framebuffer. Shared by the window and PNG export.
#pragma once

#include "ssp/app/screen.h"
#include "ssp/render/frame.h"

namespace ssp {

/// @brief Background, grid, axes, traces (first `done` columns of s.result), box, legend, title.
/// Updates s.layout, s.yshown and s.legend_hits. With `mode_label`, also writes the cmode
/// in the readout row (for static images, which have no live readout).
void compose(Framebuffer& fb, Screen& s, int done, const Theme& th, bool mode_label = false);

/// @brief Reduce synchronously and compose: the whole headless pipeline.
void render_headless(Framebuffer& fb, Screen& s, const Theme& th = {});

}  // namespace ssp
