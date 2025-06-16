#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "plot.h"

namespace py = pybind11;

PYBIND11_MODULE(signal_sniper_plot_py, m) {
    m.def("plot_buffer", [](py::array data,
                            bool is_complex,
                            bool is_float,
                            double xstart,
                            double xdelta,
                            const std::string& plot_title,
                            int line_thickness = 1,
                            std::optional<std::pair<double, double>> y_range = std::nullopt,
                            std::size_t num_traces = 1) {

        if (!data.flags() & py::array::c_style) {
            throw std::runtime_error("Input array must be C-contiguous");
        }

        py::buffer_info info = data.request();
        std::size_t num_elements = info.size;
        std::size_t elem_bytes = info.itemsize;

        xplot::plot_buffer(info.ptr,
                           num_elements,
                           elem_bytes,
                           is_complex,
                           is_float,
                           xstart,
                           xdelta,
                           plot_title,
                           line_thickness,
                           y_range,
                           num_traces);
    },
    py::arg("data"),
    py::arg("is_complex"),
    py::arg("is_float"),
    py::arg("xstart"),
    py::arg("xdelta"),
    py::arg("plot_title"),
    py::arg("line_thickness") = 1,
    py::arg("y_range") = std::nullopt,
    py::arg("num_traces") = 1,
    R"pbdoc(
        Plot an array of values (optionally complex) to an interactive X11 window.
        This wraps the xplot::plot_buffer C++ function.

        Parameters:
            data (np.ndarray): Input signal data, either interleaved complex or real
            is_complex (bool): Whether the data is complex (interleaved real/imag)
            is_float (bool): Whether the underlying data is float or int16
            xstart (float): Start time/index for x-axis
            xdelta (float): Time delta per sample
            plot_title (str): Window title
            line_thickness (int): Thickness of rendered trace
            y_range (tuple or None): (ymin, ymax) override
            num_traces (int): Number of traces in the buffer
    )pbdoc");
}
