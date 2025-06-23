#pragma once

#include <vector>
#include <cstddef>
#include <string>

namespace xplot {

enum class PlotRenderStyle { LINES, DOTS };

struct PlotSample {
    double time;
    double real;
    double imag;
};

struct Trace {
    std::vector<PlotSample> samples;
    double xstart = 0.0;
    double xdelta = 1.0;
    PlotRenderStyle style = PlotRenderStyle::LINES;
    bool visible = true;
    std::string label = "Trace";
};

Trace process_buffer(const void* data,
                     std::size_t num_elements,
                     std::size_t element_size_bytes,
                     bool is_complex,
                     bool is_float,
                     double xstart,
                     double xdelta);

}  // namespace xplot
