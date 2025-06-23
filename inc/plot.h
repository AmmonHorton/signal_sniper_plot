#pragma once

#include <string>
#include <vector>
#include <optional>
#include "plot_session.h"

namespace xplot {

// Core API for plotting memory buffers as traces
void plot_buffer(const void* data,
                 std::size_t num_elements,
                 std::size_t elem_bytes,
                 bool is_complex,
                 bool is_float,
                 double xstart,
                 double xdelta,
                 const std::string& plot_title,
                 int line_thickness = 1);

// Extended API for multiple traces with custom options
void plot_buffer_traces(const std::vector<Trace>& traces,
                        const std::string& plot_title,
                        int line_thickness = 1,
                        std::optional<std::pair<double, double>> y_range = std::nullopt,
                        std::optional<std::pair<double, double>> x_range = std::nullopt,
                        XAxisMode xmode = XAxisMode::TIME,
                        PlotMode pmode = PlotMode::Magnitude);

}  // namespace xplot
