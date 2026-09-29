// Python bindings for signal_sniper_plot (ssp).
//
//   import signal_sniper_plot as ssp
//   ssp.plot(x)                                     # any numpy dtype, real or complex
//   ssp.plot(a, b, xdelta=1/fs, names=["rx", "tx"], cmode="mag", yrange=(-1, 1))
//   ssp.plot(np.vstack([a, b]))                     # 2-D: one trace per row
//   ssp.plot(ssp.Signal(a, xstart=100, xdelta=.02, name="rx", style="dots"), b)
//   ssp.save_png("out.png", x, title="...")         # headless
//
// Arrays are read in place (no copy) and kept alive for the duration of the call.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cctype>
#include <optional>
#include <string>
#include <vector>

#include "ssp/plot.h"

namespace py = pybind11;

namespace {

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

ssp::CMode parse_cmode(const std::string& name) {
    const std::string s = lower(name);
    if (s == "auto") return ssp::CMode::Auto;
    if (s == "mag" || s == "ma" || s == "magnitude") return ssp::CMode::Mag;
    if (s == "phase" || s == "ph") return ssp::CMode::Phase;
    if (s == "real" || s == "re") return ssp::CMode::Real;
    if (s == "imag" || s == "im" || s == "imaginary") return ssp::CMode::Imag;
    if (s == "ir") return ssp::CMode::IR;
    if (s == "10log" || s == "lo" || s == "log10") return ssp::CMode::Log10;
    if (s == "20log" || s == "l2" || s == "log20") return ssp::CMode::Log20;
    throw py::value_error("unknown cmode '" + name +
                          "' (auto, mag, phase, real, imag, 10log, 20log)");
}

ssp::Style parse_style(const std::string& name) {
    const std::string s = lower(name);
    if (s == "lines") return ssp::Style::Lines;
    if (s == "dots") return ssp::Style::Dots;
    if (s == "both" || s == "linesdots") return ssp::Style::LinesDots;
    throw py::value_error("unknown style '" + name + "' (lines, dots, both)");
}

ssp::PhaseUnits parse_phunits(const std::string& name) {
    const std::string s = lower(name);
    if (s == "rad" || s == "radians" || s == "r") return ssp::PhaseUnits::Radians;
    if (s == "deg" || s == "degrees" || s == "d") return ssp::PhaseUnits::Degrees;
    if (s == "cycles" || s == "c") return ssp::PhaseUnits::Cycles;
    throw py::value_error("unknown phunits '" + name + "' (rad, deg, cycles)");
}

std::optional<ssp::Range> parse_range(const std::optional<std::pair<double, double>>& r) {
    if (!r) return std::nullopt;
    return ssp::Range{r->first, r->second};
}

std::optional<uint32_t> parse_color(const py::object& c) {
    if (c.is_none()) return std::nullopt;
    if (py::isinstance<py::int_>(c)) return c.cast<uint32_t>() & 0xFFFFFF;
    std::string s = c.cast<std::string>();
    if (!s.empty() && s[0] == '#') s.erase(0, 1);
    if (s.size() != 6) throw py::value_error("color must be 0xRRGGBB or '#RRGGBB'");
    return static_cast<uint32_t>(std::stoul(s, nullptr, 16));
}

/// A trace as the user describes it in Python; converted to ssp::Signal at call time.
struct PySignal {
    py::object data;
    double xstart = 0.0;
    double xdelta = 1.0;
    std::string name;
    std::string style = "lines";
    py::object color = py::none();
    int thickness = 0;
    bool visible = true;
};

/// Holds the arrays referenced by Signals for the duration of a call.
struct Converted {
    std::vector<ssp::Signal> signals;
    std::vector<py::object> keep;
};

ssp::DType dtype_of(const py::dtype& dt, bool& is_complex) {
    const char kind = dt.kind();
    const auto size = dt.itemsize();
    is_complex = kind == 'c';
    if (kind == 'b' || (kind == 'u' && size == 1)) return ssp::DType::U8;
    if (kind == 'i') {
        switch (size) {
            case 1: return ssp::DType::I8;
            case 2: return ssp::DType::I16;
            case 4: return ssp::DType::I32;
            case 8: return ssp::DType::I64;
        }
    }
    if (kind == 'u') {
        switch (size) {
            case 2: return ssp::DType::U16;
            case 4: return ssp::DType::U32;
            case 8: return ssp::DType::U64;
        }
    }
    if (kind == 'f' && size == 4) return ssp::DType::F32;
    if (kind == 'f' && size == 8) return ssp::DType::F64;
    if (kind == 'c' && size == 8) return ssp::DType::F32;
    if (kind == 'c' && size == 16) return ssp::DType::F64;
    throw py::type_error("unsupported dtype " + py::repr(dt).cast<std::string>());
}

/// `obj` as a 1-D or 2-D numpy array we can read in place: byte-swapped arrays and negative
/// or odd strides get one native, contiguous copy.
py::array readable_array(py::object obj) {
    py::module_ np = py::module_::import("numpy");
    py::array arr = np.attr("asarray")(obj);
    if (arr.ndim() == 0 || arr.ndim() > 2) {
        throw py::value_error("expected a 1-D or 2-D array, got " + std::to_string(arr.ndim()) + "-D");
    }
    const bool native = arr.dtype().attr("isnative").cast<bool>();
    bool strides_ok = true;
    for (py::ssize_t d = 0; d < arr.ndim(); ++d) {
        strides_ok = strides_ok && arr.strides(d) > 0 && arr.strides(d) % arr.itemsize() == 0;
    }
    if (!native || !strides_ok) {
        arr = np.attr("ascontiguousarray")(arr, py::arg("dtype") = arr.dtype().attr("newbyteorder")("="));
    }
    return arr;
}

/// numpy array (1-D or 2-D) → one Signal per row, reading memory in place when possible.
void add_array(Converted& out, py::object obj, const PySignal& spec) {
    py::array arr = readable_array(obj);
    bool cplx = false;
    const ssp::DType dt = dtype_of(arr.dtype(), cplx);
    out.keep.push_back(arr);

    const py::ssize_t rows = arr.ndim() == 2 ? arr.shape(0) : 1;
    const py::ssize_t n = arr.ndim() == 2 ? arr.shape(1) : arr.shape(0);
    const py::ssize_t elem_stride = arr.strides(arr.ndim() - 1) / arr.itemsize();
    const auto* base = static_cast<const char*>(arr.data());
    for (py::ssize_t r = 0; r < rows; ++r) {
        const char* p = base + (arr.ndim() == 2 ? r * arr.strides(0) : 0);
        ssp::Signal s(p, dt, cplx, static_cast<std::size_t>(n), spec.xstart, spec.xdelta);
        s.stride = elem_stride;
        s.name = spec.name.empty() || rows == 1 ? spec.name : spec.name + "[" + std::to_string(r) + "]";
        s.style = parse_style(spec.style);
        s.color = parse_color(spec.color);
        s.thickness = spec.thickness;
        s.visible = spec.visible;
        out.signals.push_back(std::move(s));
    }
}

/// Positional args: arrays, Signals, or lists/tuples of those.
Converted convert(const py::args& data, double xstart, double xdelta,
                  const std::optional<std::vector<std::string>>& names) {
    Converted out;
    std::vector<py::object> items;
    for (const py::handle h : data) {
        if (py::isinstance<py::list>(h) || py::isinstance<py::tuple>(h)) {
            for (const py::handle e : h) items.push_back(py::reinterpret_borrow<py::object>(e));
        } else {
            items.push_back(py::reinterpret_borrow<py::object>(h));
        }
    }
    if (items.empty()) throw py::value_error("nothing to plot: pass at least one array");
    for (const auto& item : items) {
        if (py::isinstance<PySignal>(item)) {
            add_array(out, item.cast<const PySignal&>().data, item.cast<const PySignal&>());
        } else {
            PySignal spec;
            spec.xstart = xstart;
            spec.xdelta = xdelta;
            add_array(out, item, spec);
        }
    }
    if (names) {
        if (names->size() != out.signals.size()) {
            throw py::value_error("names has " + std::to_string(names->size()) + " entries for " +
                                  std::to_string(out.signals.size()) + " traces");
        }
        for (std::size_t i = 0; i < names->size(); ++i) out.signals[i].name = (*names)[i];
    }
    return out;
}

ssp::PlotOptions make_options(const std::string& title, const std::string& cmode,
                              const std::optional<std::pair<double, double>>& xrange,
                              const std::optional<std::pair<double, double>>& yrange, bool index,
                              int thickness, bool grid, bool legend, const std::string& phunits,
                              int width, int height) {
    ssp::PlotOptions o;
    o.title = title;
    o.cmode = parse_cmode(cmode);
    o.xrange = parse_range(xrange);
    o.yrange = parse_range(yrange);
    o.index = index;
    o.thickness = thickness;
    o.grid = grid;
    o.legend = legend;
    o.phunits = parse_phunits(phunits);
    o.width = width;
    o.height = height;
    return o;
}

/// Blocking window with the GIL released; Ctrl-C closes the window and raises KeyboardInterrupt.
void run_window(const Converted& c, const ssp::PlotOptions& o) {
    ssp::Session session(c.signals, o);
    bool interrupted = false;
    session.set_interrupt_check([&] {
        py::gil_scoped_acquire gil;
        if (PyErr_CheckSignals() != 0) interrupted = true;
        return interrupted;
    });
    {
        py::gil_scoped_release nogil;
        session.run();
    }
    if (interrupted) throw py::error_already_set();
}

ssp::Colormap parse_cmap(const std::string& name) {
    std::string s = lower(name);
    s.erase(std::remove_if(s.begin(), s.end(), [](char c) { return c == ' ' || c == '_' || c == '-'; }), s.end());
    if (s == "greyscale" || s == "grayscale" || s == "grey" || s == "gray") return ssp::Colormap::Greyscale;
    if (s == "ramp") return ssp::Colormap::Ramp;
    if (s == "colorwheel" || s == "wheel") return ssp::Colormap::ColorWheel;
    if (s == "spectrum") return ssp::Colormap::Spectrum;
    if (s == "calewhite") return ssp::Colormap::CalEWhite;
    if (s == "hotdesat") return ssp::Colormap::HotDesat;
    if (s == "sunset") return ssp::Colormap::Sunset;
    if (s == "hot") return ssp::Colormap::Hot;
    if (s == "cold") return ssp::Colormap::Cold;
    throw py::value_error("unknown cmap '" + name +
                          "' (greyscale, ramp, colorwheel, spectrum, calewhite, hotdesat, sunset, hot, cold)");
}

ssp::Reduce parse_reduce(const std::string& name) {
    const std::string s = lower(name);
    if (s == "max") return ssp::Reduce::Max;
    if (s == "min") return ssp::Reduce::Min;
    if (s == "mean" || s == "avg") return ssp::Reduce::Mean;
    if (s == "maxabs") return ssp::Reduce::MaxAbs;
    if (s == "first") return ssp::Reduce::First;
    throw py::value_error("unknown reduce '" + name + "' (max, min, mean, maxabs, first)");
}

/// A raster: 2-D array (rows = frames) or 1-D array cut into frames of `subsize`.
struct RasterInput {
    ssp::Signal signal;
    ssp::RasterOptions options;
    py::array keep;
};

RasterInput make_raster(py::object data, std::optional<std::size_t> subsize, double xstart,
                        double xdelta, double ystart, double ydelta, const std::string& title,
                        const std::string& cmode, const std::optional<std::pair<double, double>>& xrange,
                        const std::optional<std::pair<double, double>>& yrange,
                        const std::optional<std::pair<double, double>>& zrange, const std::string& cmap,
                        const std::string& reduce, bool index, bool grid, const std::string& phunits,
                        int width, int height) {
    RasterInput in;
    in.keep = readable_array(std::move(data));
    const py::array& a = in.keep;
    bool cplx = false;
    const ssp::DType dt = dtype_of(a.dtype(), cplx);
    ssp::RasterOptions& o = in.options;
    std::size_t n = 0;
    if (a.ndim() == 2) {
        if (subsize && *subsize != static_cast<std::size_t>(a.shape(1))) {
            throw py::value_error("subsize must match the array's row length (or be omitted)");
        }
        o.subsize = static_cast<std::size_t>(a.shape(1));
        o.frame_stride = a.strides(0) / a.itemsize();
        n = static_cast<std::size_t>(a.shape(0)) * o.subsize;
    } else {
        if (!subsize) throw py::value_error("a 1-D array needs subsize= (samples per frame)");
        o.subsize = *subsize;
        n = static_cast<std::size_t>(a.shape(0));
    }
    in.signal = ssp::Signal(a.data(), dt, cplx, n, xstart, xdelta);
    in.signal.stride = a.strides(a.ndim() - 1) / a.itemsize();
    o.title = title;
    o.cmode = parse_cmode(cmode);
    o.phunits = parse_phunits(phunits);
    o.ystart = ystart;
    o.ydelta = ydelta;
    o.xrange = parse_range(xrange);
    o.yrange = parse_range(yrange);
    o.zrange = parse_range(zrange);
    o.cmap = parse_cmap(cmap);
    o.reduce = parse_reduce(reduce);
    o.index = index;
    o.grid = grid;
    o.width = width;
    o.height = height;
    return in;
}

#define SSP_RASTER_KWARGS                                                                      \
    py::kw_only(), py::arg("subsize") = py::none(), py::arg("xstart") = 0.0,                   \
        py::arg("xdelta") = 1.0, py::arg("ystart") = 0.0, py::arg("ydelta") = 1.0,              \
        py::arg("title") = "", py::arg("cmode") = "auto", py::arg("xrange") = py::none(),        \
        py::arg("yrange") = py::none(), py::arg("zrange") = py::none(), py::arg("cmap") = "ramp", \
        py::arg("reduce") = "max", py::arg("index") = false, py::arg("grid") = false,           \
        py::arg("phunits") = "rad", py::arg("width") = 1000, py::arg("height") = 600

#define SSP_PLOT_KWARGS                                                                    \
    py::kw_only(), py::arg("title") = "", py::arg("cmode") = "auto", py::arg("xstart") = 0.0, \
        py::arg("xdelta") = 1.0, py::arg("names") = py::none(), py::arg("xrange") = py::none(), \
        py::arg("yrange") = py::none(), py::arg("index") = false, py::arg("thickness") = 1,   \
        py::arg("grid") = true, py::arg("legend") = true, py::arg("phunits") = "rad"

}  // namespace

PYBIND11_MODULE(signal_sniper_plot, m) {
    m.doc() = "Fast interactive DSP plots (xplot) for numpy arrays of any size.";

    py::class_<PySignal>(m, "Signal", R"doc(
        One trace with its own x axis and style.

        Signal(data, xstart=0.0, xdelta=1.0, name="", style="lines", color=None,
               thickness=0, visible=True)

        style: "lines", "dots" or "both". color: 0xRRGGBB or "#RRGGBB".
        thickness: pixels, 0 = the plot's default.
    )doc")
        .def(py::init([](py::object data, double xstart, double xdelta, std::string name,
                         std::string style, py::object color, int thickness, bool visible) {
                 parse_style(style);  // validate early
                 parse_color(color);
                 return PySignal{std::move(data), xstart, xdelta, std::move(name),
                                 std::move(style), std::move(color), thickness, visible};
             }),
             py::arg("data"), py::arg("xstart") = 0.0, py::arg("xdelta") = 1.0,
             py::arg("name") = "", py::arg("style") = "lines", py::arg("color") = py::none(),
             py::arg("thickness") = 0, py::arg("visible") = true)
        .def_readwrite("data", &PySignal::data)
        .def_readwrite("xstart", &PySignal::xstart)
        .def_readwrite("xdelta", &PySignal::xdelta)
        .def_readwrite("name", &PySignal::name)
        .def_readwrite("style", &PySignal::style)
        .def_readwrite("color", &PySignal::color)
        .def_readwrite("thickness", &PySignal::thickness)
        .def_readwrite("visible", &PySignal::visible);

    m.def(
        "plot",
        [](const py::args& data, const std::string& title, const std::string& cmode,
           double xstart, double xdelta, const std::optional<std::vector<std::string>>& names,
           const std::optional<std::pair<double, double>>& xrange,
           const std::optional<std::pair<double, double>>& yrange, bool index, int thickness,
           bool grid, bool legend, const std::string& phunits, int width, int height) {
            const Converted c = convert(data, xstart, xdelta, names);
            run_window(c, make_options(title, cmode, xrange, yrange, index, thickness, grid,
                                       legend, phunits, width, height));
        },
        SSP_PLOT_KWARGS, py::arg("width") = 1000, py::arg("height") = 600,
        R"doc(
        Open an interactive plot window; returns when it is closed.

        Positional arguments are numpy arrays (1-D: one trace, 2-D: one trace per row),
        Signal objects, or lists of those. xstart/xdelta apply to bare arrays.

        cmode: auto, mag, phase, real, imag, 10log, 20log.  xrange/yrange: (lo, hi).
        In the window: ? shows every key, middle-click (or m) opens the menu. Drag to zoom,
        right-click to unzoom, click the legend to show/hide (right-click: lines/dots/both).
    )doc");

    m.def(
        "save_png",
        [](const std::string& path, const py::args& data, const std::string& title,
           const std::string& cmode, double xstart, double xdelta,
           const std::optional<std::vector<std::string>>& names,
           const std::optional<std::pair<double, double>>& xrange,
           const std::optional<std::pair<double, double>>& yrange, bool index, int thickness,
           bool grid, bool legend, const std::string& phunits, int width, int height) {
            const Converted c = convert(data, xstart, xdelta, names);
            const ssp::PlotOptions o = make_options(title, cmode, xrange, yrange, index, thickness,
                                                    grid, legend, phunits, width, height);
            py::gil_scoped_release nogil;
            ssp::save_png(c.signals, o, path, width, height);
        },
        py::arg("path"), SSP_PLOT_KWARGS, py::arg("width") = 1000, py::arg("height") = 600,
        "Render the plot without a display and write it to `path` as PNG. Same arguments as plot().");

    m.def(
        "raster",
        [](py::object data, std::optional<std::size_t> subsize, double xstart, double xdelta,
           double ystart, double ydelta, const std::string& title, const std::string& cmode,
           const std::optional<std::pair<double, double>>& xrange,
           const std::optional<std::pair<double, double>>& yrange,
           const std::optional<std::pair<double, double>>& zrange, const std::string& cmap,
           const std::string& reduce, bool index, bool grid, const std::string& phunits,
           int width, int height) {
            const RasterInput in = make_raster(std::move(data), subsize, xstart, xdelta, ystart, ydelta,
                                               title, cmode, xrange, yrange, zrange, cmap, reduce,
                                               index, grid, phunits, width, height);
            ssp::Session session(in.signal, in.options);
            bool interrupted = false;
            session.set_interrupt_check([&] {
                py::gil_scoped_acquire gil;
                if (PyErr_CheckSignals() != 0) interrupted = true;
                return interrupted;
            });
            {
                py::gil_scoped_release nogil;
                session.run();
            }
            if (interrupted) throw py::error_already_set();
        },
        py::arg("data"), SSP_RASTER_KWARGS,
        R"doc(
        Open an interactive raster (xraster) window; returns when it is closed.

        data: 2-D array (one frame per row, frame 0 drawn at the top), or a 1-D array with
        subsize= samples per frame. Columns use xstart/xdelta, frames ystart/ydelta.
        cmode: auto, mag, phase, real, imag, 10log, 20log.  zrange: colour range (lo, hi).
        cmap: greyscale, ramp, colorwheel, spectrum, calewhite, hotdesat, sunset, hot, cold.
        reduce: how a pixel combines the samples under it: max, min, mean, maxabs, first.
        In the window: x / y cut the row / column under the pointer (Esc returns),
        c changes colormap, [ ] slide the colour range, ? shows every key.
    )doc");

    m.def(
        "save_raster_png",
        [](const std::string& path, py::object data, std::optional<std::size_t> subsize,
           double xstart, double xdelta, double ystart, double ydelta, const std::string& title,
           const std::string& cmode, const std::optional<std::pair<double, double>>& xrange,
           const std::optional<std::pair<double, double>>& yrange,
           const std::optional<std::pair<double, double>>& zrange, const std::string& cmap,
           const std::string& reduce, bool index, bool grid, const std::string& phunits,
           int width, int height) {
            const RasterInput in = make_raster(std::move(data), subsize, xstart, xdelta, ystart, ydelta,
                                               title, cmode, xrange, yrange, zrange, cmap, reduce,
                                               index, grid, phunits, width, height);
            py::gil_scoped_release nogil;
            ssp::save_raster_png(in.signal, in.options, path, width, height);
        },
        py::arg("path"), py::arg("data"), SSP_RASTER_KWARGS,
        "Render a raster without a display and write it to `path` as PNG. Same arguments as raster().");

    // Compatibility with the original API; prefer plot().
    m.def(
        "plot_buffer",
        [](py::object data, double xstart, double xdelta, const std::string& plot_title,
           int line_thickness, const std::optional<std::pair<double, double>>& y_range,
           const std::optional<std::pair<double, double>>& x_range, std::size_t num_traces) {
            if (num_traces == 0) throw py::value_error("num_traces must be > 0");
            py::array arr = py::module_::import("numpy").attr("asarray")(data);
            if (num_traces > 1) {
                const auto total = static_cast<std::size_t>(arr.size());
                if (total % num_traces != 0) {
                    throw py::value_error("array length is not a multiple of num_traces");
                }
                arr = arr.attr("reshape")(num_traces, total / num_traces);
            }
            const Converted c = convert(py::make_tuple(arr), xstart, xdelta, std::nullopt);
            run_window(c, make_options(plot_title, "auto", x_range, y_range, false, line_thickness,
                                       true, true, "rad", 1000, 600));
        },
        py::arg("data"), py::arg("xstart") = 0.0, py::arg("xdelta") = 1.0,
        py::arg("plot_title") = "Plot", py::arg("line_thickness") = 2,
        py::arg("y_range") = py::none(), py::arg("x_range") = py::none(),
        py::arg("num_traces") = 1,
        "Deprecated: use plot(). num_traces splits the array into that many equal, "
        "consecutive traces.");
}
