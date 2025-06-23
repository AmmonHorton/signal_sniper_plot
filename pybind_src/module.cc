#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "plot.h"
#include "trace_utils.h"

namespace py = pybind11;

bool is_numpy_complex(const py::dtype& dt) {
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
                            std::optional<std::pair<double, double>> x_range,
                            std::size_t num_traces) {
        
        if (!(data.flags() & py::array::c_style)) {
            throw std::runtime_error("Input array must be C-contiguous");
        }

        py::buffer_info info = data.request();
        std::size_t total_elements = info.size;
        std::size_t elem_bytes = info.itemsize;

        if (num_traces == 0) throw std::runtime_error("num_traces must be > 0");

        py::dtype dtype = data.dtype();
        bool is_complex = is_numpy_complex(dtype);
        bool is_float = is_numpy_float(dtype) || is_complex;

        std::size_t trace_len = total_elements / num_traces;
        std::vector<xplot::Trace> traces;

        const char* base_ptr = static_cast<const char*>(info.ptr);

        for (std::size_t t = 0; t < num_traces; ++t) {
            const void* trace_ptr = base_ptr + t * trace_len * elem_bytes;
            xplot::Trace trace = xplot::process_buffer(
                trace_ptr, trace_len, elem_bytes,
                is_complex, is_float, xstart, xdelta);
            trace.label = "Trace " + std::to_string(t);
            traces.push_back(std::move(trace));
        }

        xplot::plot_buffer_traces(
            traces,
            plot_title,
            line_thickness,
            y_range,
            x_range,
            xplot::XAxisMode::TIME,
            xplot::PlotMode::Magnitude
        );
    },
    py::arg("data"),
    py::arg("xstart") = 0.0,
    py::arg("xdelta") = 1.0,
    py::arg("plot_title") = "Plot",
    py::arg("line_thickness") = 2,
    py::arg("y_range") = std::nullopt,
    py::arg("x_range") = std::nullopt,
    py::arg("num_traces") = 1,
    R"pbdoc(
        Plot one or more interleaved signals from a 1D NumPy array to an X11 window.

        Parameters:
            data (np.ndarray): Real or complex 1D array (float32, int16, complex64, etc).
            xstart (float): Start value of x-axis.
            xdelta (float): Step between x-axis samples.
            plot_title (str): Title for the plot window.
            line_thickness (int): Line or dot size.
            y_range (Optional[Tuple[float, float]]): Fixed y-limits.
            num_traces (int): Number of interleaved traces.
    )pbdoc");
}
