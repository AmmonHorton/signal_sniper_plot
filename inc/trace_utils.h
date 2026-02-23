/// @file trace_utils.h
/// @brief Core sample and trace types, and the raw-buffer conversion utility.
#pragma once

#include <vector>
#include <cstddef>
#include <string>

namespace xplot {

/// @brief Rendering style for a trace.
enum class PlotRenderStyle {
    LINES, ///< Connect samples with continuous line segments.
    DOTS   ///< Draw each sample as an individual filled circle.
};

/// @brief A single time-domain sample with real and imaginary components.
struct PlotSample {
    double time; ///< X-axis coordinate (time or index-derived value).
    double real; ///< Real component of the sample value.
    double imag; ///< Imaginary component; zero for real-only signals.
};

/// @brief A named, styled collection of PlotSamples rendered as one trace.
struct Trace {
    std::vector<PlotSample> samples;                ///< Ordered sample data.
    double xstart = 0.0;                            ///< X value of the first sample.
    double xdelta = 1.0;                            ///< X-axis step between adjacent samples.
    PlotRenderStyle style = PlotRenderStyle::LINES; ///< How the trace is drawn.
    bool visible = true;                            ///< Whether the trace is currently shown.
    std::string label = "Trace";                    ///< Display name shown in the legend.
};

/// @brief Convert a raw memory buffer into a Trace of PlotSamples.
/// @param data               Pointer to the raw sample data.
/// @param num_elements       Number of elements in the buffer.
/// @param element_size_bytes Size of each element in bytes.
/// @param is_complex         True if elements are interleaved real/imag pairs.
/// @param is_float           True if elements are floating-point; false for integer.
/// @param xstart             X-axis value assigned to the first sample.
/// @param xdelta             X-axis increment between consecutive samples.
/// @return A Trace with @p num_elements PlotSamples populated from the buffer.
Trace process_buffer(const void* data,
                     std::size_t num_elements,
                     std::size_t element_size_bytes,
                     bool is_complex,
                     bool is_float,
                     double xstart,
                     double xdelta);

}  // namespace xplot
