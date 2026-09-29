/// @file plot.h
/// @brief Public entry points for xplot.
///
/// Three layers, each only filling in defaults for the one below:
///   ssp::plot(signal)                       one trace, all defaults
///   ssp::plot(signals, options)             any number of traces, every option
///   ssp::Session                            full control (interrupt hook; later: live data)
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "ssp/types.h"

namespace ssp {

/// @brief Plot-wide options. Every field is optional.
struct PlotOptions {
    std::string title;
    CMode cmode = CMode::Auto;
    PhaseUnits phunits = PhaseUnits::Radians;
    std::optional<Range> xrange;  ///< Fixed x view; autoscaled from the data when empty.
    std::optional<Range> yrange;  ///< Fixed y view; autoscaled from the data in view when empty.
    bool index = false;           ///< Use the sample index (0-based) as x instead of xstart/xdelta.
    int thickness = 1;            ///< Default line thickness in pixels.
    bool grid = true;
    bool legend = true;
    int width = 1000;             ///< Initial window size in pixels.
    int height = 600;
};

/// @brief xraster options. `subsize` is required; everything else is optional.
///
/// The data is a Signal of frames: cell (row r, column c) is the sample at offset
/// r * frame_stride + c * signal.stride from signal.data (offsets count samples; a complex pair
/// is one sample), and exists while r * subsize + c < signal.n, so a partial last frame shows
/// as empty cells.
/// Columns take the signal's xstart/xdelta; rows (frames) take ystart/ydelta.
struct RasterOptions {
    std::string title;
    CMode cmode = CMode::Auto;       ///< Anything but IR; Auto = Real, or Mag for complex data.
    PhaseUnits phunits = PhaseUnits::Radians;
    std::size_t subsize = 0;         ///< Samples per frame (row).
    std::ptrdiff_t frame_stride = 0; ///< Samples between frame starts; 0 = subsize * stride.
    double ystart = 0.0;             ///< y of frame 0 (drawn at the top).
    double ydelta = 1.0;             ///< y step between frames; must be > 0.
    std::optional<Range> xrange, yrange;
    std::optional<Range> zrange;     ///< Colour range; autoscaled from the data in view when empty.
    Colormap cmap = Colormap::Ramp;
    Reduce reduce = Reduce::Max;     ///< How a pixel combines the samples under it.
    bool index = false;              ///< Axes in column/row numbers.
    bool grid = false;
    int width = 1000;
    int height = 600;
};

/// @brief An interactive plot window.
class Session {
public:
    /// @throws std::invalid_argument for unusable signals or options.
    explicit Session(std::vector<Signal> signals, PlotOptions options = {});
    /// @brief An xraster window.
    /// @throws std::invalid_argument for unusable data or options.
    Session(Signal data, RasterOptions options);

    /// @brief Called about every 100 ms while the window is open; returning true closes it
    /// (Python uses this for Ctrl-C).
    void set_interrupt_check(std::function<bool()> fn) { interrupt_ = std::move(fn); }

    /// @brief Open the window and block until it is closed.
    /// @throws std::runtime_error when no display is available.
    void run();

private:
    std::vector<Signal> signals_;
    PlotOptions options_;
    std::optional<RasterOptions> raster_;
    std::function<bool()> interrupt_;
};

/// @brief Open an interactive plot window and block until it is closed.
void plot(const std::vector<Signal>& signals, const PlotOptions& options = {});
inline void plot(const Signal& signal, const PlotOptions& options = {}) {
    plot(std::vector<Signal>{signal}, options);
}

/// @brief Open an interactive raster window and block until it is closed.
/// Keys x/y open a cut of the row/column under the pointer; Esc returns.
void raster(const Signal& data, const RasterOptions& options);

/// @brief Render a raster headlessly and write it as a PNG.
void save_raster_png(const Signal& data, const RasterOptions& options, const std::string& path,
                     int width = 1000, int height = 600);

/// @brief Render the plot headlessly (no display needed) and write it as a PNG.
/// @throws std::invalid_argument for unusable input, std::runtime_error on I/O failure.
void save_png(const std::vector<Signal>& signals, const PlotOptions& options,
              const std::string& path, int width = 1000, int height = 600);

}  // namespace ssp
