/// @file plot.h
/// @brief Public API for opening an interactive plot window from raw buffers or pre-built traces.
#pragma once

#include <string>
#include <vector>
#include <optional>
#include "plot_session.h"

namespace xplot {

/// @brief Plot a single raw buffer as an interactive trace window.
/// Convenience wrapper that calls process_buffer() and opens a PlotSession.
/// @param data           Pointer to the raw sample data.
/// @param num_elements   Number of elements in the buffer.
/// @param elem_bytes     Size of each element in bytes.
/// @param is_complex     True if elements are interleaved real/imag pairs.
/// @param is_float       True if elements are floating-point; false for integer.
/// @param xstart         X-axis start value.
/// @param xdelta         X-axis step between samples.
/// @param plot_title     Window title string.
/// @param line_thickness Pixel thickness of the plotted line.
void plot_buffer(const void* data,
                 std::size_t num_elements,
                 std::size_t elem_bytes,
                 bool is_complex,
                 bool is_float,
                 double xstart,
                 double xdelta,
                 const std::string& plot_title,
                 int line_thickness = 1);

/// @brief Plot one or more pre-built traces in an interactive window.
/// Blocks until the window is closed by the user.
/// @param traces         Ordered list of traces to display.
/// @param plot_title     Window title string.
/// @param line_thickness Pixel thickness of plotted lines.
/// @param y_range        Optional fixed Y-axis bounds; auto-scales when empty.
/// @param x_range        Optional fixed X-axis bounds; auto-scales when empty.
/// @param xmode          X-axis label mode (TIME or INDEX).
/// @param pmode          Initial plot mode (Magnitude, Real, Imag, Phase, IQ).
void plot_buffer_traces(const std::vector<Trace>& traces,
                        const std::string& plot_title,
                        int line_thickness = 1,
                        std::optional<std::pair<double, double>> y_range = std::nullopt,
                        std::optional<std::pair<double, double>> x_range = std::nullopt,
                        XAxisMode xmode = XAxisMode::TIME,
                        PlotMode pmode = PlotMode::Magnitude);

}  // namespace xplot
