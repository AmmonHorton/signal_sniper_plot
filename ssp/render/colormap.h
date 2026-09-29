/// @file colormap.h
/// @brief 256-entry colour lookup tables for rasters.
#pragma once

#include <array>
#include <cstdint>

#include "ssp/types.h"

namespace ssp {

using ColorLut = std::array<uint32_t, 256>;  ///< 0x00RRGGBB, low z first.

const ColorLut& colormap_lut(Colormap c);
const char* colormap_name(Colormap c);

/// @brief LUT index for z in [zmin, zmax] (clamped); -1 for NaN.
inline int lut_index(double z, double zmin, double zmax) {
    if (!(z == z)) return -1;
    const double f = (z - zmin) / (zmax - zmin) * 256.0;
    return f <= 0.0 ? 0 : f >= 255.0 ? 255 : static_cast<int>(f);
}

}  // namespace ssp
