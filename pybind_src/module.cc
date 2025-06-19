#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <iostream>
#include "plot.h"

namespace py = pybind11;

bool is_numpy_complex(const py::dtype& dt) {
    std::cout << dt.kind() << "\n";
    return dt.kind() == 'c';  // complex64 or complex128
}

bool is_numpy_float(const py::dtype& dt) {
    return dt.kind() == 'f';  // float32 or float64
}

PYBIND11_MODULE(signal_sniper_plot_py, m) {
    m.doc() = "Python bindings for signal_sniper_plot using Pybind11";

    m.def("plot_buffer", [](py::array data,
                            double xstart,
                            double xdelta,
                            const std::string& plot_title,
                            int line_thickness,
                            std::optional<std::pair<double, double>> y_range,
                            std::size_t num_traces) {
        
        if (!(data.flags() & py::array::c_style)) {
            throw std::runtime_error("Input array must be C-contiguous");
        }

        py::buffer_info info = data.request();
        std::size_t num_elements = info.size;
        std::size_t elem_bytes = info.itemsize;

        py::dtype dtype = data.dtype();
        bool is_complex = is_numpy_complex(dtype);
        bool is_float = is_numpy_float(dtype) || is_complex;

        if (is_complex) {
            elem_bytes /= 2;
        }

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
    py::arg("xstart") = 0.0,
    py::arg("xdelta") = 1.0,
    py::arg("plot_title") = "Plot",
    py::arg("line_thickness") = 2,
    py::arg("y_range") = std::nullopt,
    py::arg("num_traces") = 1,
    R"pbdoc(
        Plot a 1D real or complex NumPy array to an interactive X11 window.

        Parameters:
            data (np.ndarray): Input signal data. Can be float32, int16, complex64, etc.
            xstart (float): X-axis start value. Default = 0.0
            xdelta (float): X-axis increment per sample. Default = 1.0
            plot_title (str): Window title. Default = "Plot"
            line_thickness (int): Line thickness. Default = 2
            y_range (tuple[float, float] or None): Fixed y-axis limits. Default = None
            num_traces (int): Number of interleaved traces. Default = 1
    )pbdoc");
}
