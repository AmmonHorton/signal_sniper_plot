/// @file cut.h
/// @brief x- and y-cuts: an xplot screen of one raster row or column, over the raster's
/// zoom box. The cut owns its own settings and zoom stack, so nothing done in it can change
/// the raster underneath.
#pragma once

#include <memory>

#include "app/screen.h"

namespace ssp {

/// @brief A screen plotting row (xcut) or column `index` of `raster`'s data inside its
/// current zoom box.
std::unique_ptr<Screen> make_cut(const Screen& raster, bool xcut, std::size_t index);

/// @brief Move `cut` to the next/previous row or column within the box, keeping its zoom.
/// @return false (with a message) at the edge of the box.
bool step_cut(Screen& cut, int step);

}  // namespace ssp
