/// @file overlay.h
/// @brief Things drawn over the composed plot every frame: crosshair, zoom box, readout, status.
#pragma once

#include <vector>

#include "ssp/app/screen.h"
#include "ssp/render/frame.h"

namespace ssp {

/// @brief Draw overlays onto `fb` (a copy of the composed plot). Returns every rect touched so
/// the caller can restore those areas from the plot next frame and upload only them.
/// @param progress fraction of the current reduce finished, in [0, 1].
std::vector<Rect> draw_overlay(Framebuffer& fb, const Screen& s, const Theme& th, bool rendering,
                               double progress);

}  // namespace ssp
