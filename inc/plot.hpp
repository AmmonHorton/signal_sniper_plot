#pragma once
#include <cstddef>
#include <string>

namespace xplot {

void plot_buffer(
    const void* data,
    std::size_t num_elements,
    std::size_t element_size_bytes,
    bool is_complex,
    bool is_float,
    double xstart_j1950_sec,
    double xdelta_sec,
    const std::string& plot_title
);

}  // namespace xplot