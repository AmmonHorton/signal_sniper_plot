#pragma once
#include <cstddef>
#include <string>
#include <optional>

namespace xplot {

void plot_buffer(
    const void* data,
    std::size_t num_elements,
    std::size_t element_size_bytes,
    bool is_complex,
    bool is_float,
    double xstart_j1950_sec = 0.0,
    double xdelta_sec = 1.0,
    const std::string& plot_title = "Plot",
    int line_thickness = 2,
    std::optional<std::pair<double, double>> y_range = std::nullopt
);

}  // namespace xplot