/// @file plot.h
/// @brief Public entry points for xplot.
#pragma once

#include <string>
#include <vector>

#include "ssp/types.h"

namespace ssp {

/// @brief Plot-wide options. Every field is optional.
struct PlotOptions {
    std::string title;
    CMode cmode = CMode::Auto;
    PhaseUnits phunits = PhaseUnits::Radians;
    std::optional<Range> xrange;  ///< Fixed x view; autoscaled from the data when empty.
    std::optional<Range> yrange;  ///< Fixed y view; autoscaled from the data in view when empty.
    bool index = false;           ///< Use the sample index (0-based) as x instead of xstart/xdelta.
    int thickness = 1;            ///< Default line thickness in pixels.
    bool grid = true;
    bool legend = true;
};

/// @brief Render the plot headlessly (no display needed) and write it as a PNG.
/// @throws std::invalid_argument for unusable input, std::runtime_error on I/O failure.
void save_png(const std::vector<Signal>& signals, const PlotOptions& options,
              const std::string& path, int width = 1000, int height = 600);

}  // namespace ssp
