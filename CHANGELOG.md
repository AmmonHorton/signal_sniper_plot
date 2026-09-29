# Changelog

## 2.0.0

A rewrite for large data, and the start of xraster.

### Added
- **xraster:** `raster()` / `ssp::raster` for 2-D data with per-pixel reduction, XMidas
  colormaps, a colour bar and z readout, and x/y cuts (an xplot of the row/column under the
  pointer within the zoom box; `Esc` returns, `PgUp`/`PgDn` step).
- **Large data:** samples are read in place in any dtype (including `np.memmap` and strided
  views); a min/max level-of-detail pyramid makes views after the first cost O(pixels).
- **Rendering you can stop:** a worker thread renders progressively; any click stops it,
  `Space` resumes. Unzoom is instant (each zoom level's result is kept).
- **xplot:** Imag-vs-Real density mode, 10/20·log10 modes, phase units, marker with dx/dy
  readout, wheel zoom, pan, typed x/y ranges, exact dots mode, middle-click menu, `?` help.
- **Headless:** `save_png` / `save_raster_png`.
- **Packages:** manylinux_2_28 wheel (RHEL/CentOS Stream 8+, Ubuntu 20.04+), and RPM
  (RHEL 8+) / .deb (Ubuntu 20.04+) packages of the C++ library with a pkg-config file.

### Changed
- The Python module is now `signal_sniper_plot` (`import signal_sniper_plot as ssp`); the pip
  name is still `signal_sniper_plot_py`. `import signal_sniper_plot_py` still works with a
  `DeprecationWarning` until 3.0.
- C++ API: `ssp::plot`, `ssp::raster`, `ssp::Session` in `ssp/plot.h` (namespace `ssp`).

### Removed
- The 1.x C++ API (`xplot::plot_buffer`, `plot_buffer_traces`, `PlotSession`). Python's
  `plot_buffer` still works.

### Fixed (from 1.x)
- Phase used `atan2(real, imag)`; autoscale looked only at the first trace; the idle loop
  busy-polled; no X display called `exit(1)`, killing Python; `uint16` data was read as
  signed; the GIL was held while a window was open.
