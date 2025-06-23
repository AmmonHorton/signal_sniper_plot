// plot.cc
#include "plot.h"
#include "plot_session.h"
#include "trace_utils.h"
#include <iostream>

namespace xplot {

void plot_buffer(const void* data,
                 std::size_t num_elements,
                 std::size_t elem_bytes,
                 bool is_complex,
                 bool is_float,
                 double xstart,
                 double xdelta,
                 const std::string& plot_title,
                 int line_thickness) {
    if (num_elements == 0) return;

    Trace trace = process_buffer(data, num_elements, elem_bytes,
                                 is_complex, is_float, xstart, xdelta);

    trace.label = "Trace 0";

    PlotSession session(plot_title, line_thickness);
    session.add_trace(trace);
    session.run();
}

void plot_buffer_traces(const std::vector<Trace>& traces,
                        const std::string& plot_title,
                        int line_thickness,
                        std::optional<std::pair<double, double>> y_range,
                        std::optional<std::pair<double, double>> x_range,
                        XAxisMode xmode,
                        PlotMode pmode) {
    if (traces.empty()) return;

    PlotSession session(plot_title, line_thickness);

    for (const auto& trace : traces) {
        session.add_trace(trace);
    }

    session.set_axis_mode(xmode);
    session.set_plot_mode(pmode);

    if (y_range) {
        session.set_fixed_y_range(y_range);
    }

    if (x_range) {
        session.set_fixed_x_range(x_range);
    }

    session.run();
}

} // namespace xplot
