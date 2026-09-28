# Signal Sniper Plot

Fast interactive DSP plotting for C++ and Python, in the spirit of XMidas `xplot` / `xraster`
and SigPlot: line plots of 1-D signals and raster (image) plots of frames, in an X11 window,
for arrays of any size.

- **Big data.** Samples are read in place, with no copies or conversion; `np.memmap` works.
  A min/max level-of-detail pyramid makes every view after the first cost O(pixels): about
  7 ms for a 1.6 GB signal.
- **Watch it render, stop it any time.** Rendering runs on a worker thread and fills the
  window left to right (rasters: top to bottom). Any click stops it, and the same drag zooms.
- **xplot.** Modes Magnitude, Phase, Real, Imaginary, Imag-vs-Real (a density plot of a time
  window), 10·log10, 20·log10. Box, wheel and typed-range zoom; instant unzoom; pan; marker
  with dx/dy readout; lines / dots (dots light only pixels that hold samples); legend; grid.
- **xraster.** Per-pixel max / min / mean / max-abs / first reduction, XMidas colormaps,
  colour bar, and x/y cuts: `x` / `y` open a line plot of the row / column under the pointer
  within the zoom box, isolated from the raster; `Esc` returns.
- **Headless.** `save_png` / `save_raster_png` render without a display.

In the window, `?` lists every key and middle-click opens the menu.

---

## Using it

Python (see [README_pypi.md](README_pypi.md) for the full reference):

```python
import numpy as np
import signal_sniper_plot_py as ssp

ssp.plot(iq, xdelta=1 / fs, cmode="mag")
ssp.plot(ssp.Signal(rx, xdelta=1 / fs, name="rx"), ssp.Signal(tx, xstart=2e-3, xdelta=1 / fs, name="tx"))
ssp.raster(spectrogram, cmode="20log", cmap="ramp")
```

C++ (`#include "ssp/plot.h"`, Bazel target `//:ssp`):

```cpp
std::vector<std::complex<float>> iq = ...;
ssp::PlotOptions o;                 // every field is optional
o.title = "capture";
o.yrange = ssp::Range{-1.5, 1.5};
ssp::plot({ssp::Signal(iq, 0.0, 1e-7, "iq")}, o);   // blocks until the window closes

std::vector<float> spec = ...;      // frames laid end to end
ssp::RasterOptions r;
r.subsize = 4096;                   // samples per frame
r.cmode = ssp::CMode::Log20;
ssp::raster(ssp::Signal(spec), r);
```

`ssp::Session` is the full-control layer under `plot` and `raster` (e.g. an interrupt hook).
The 1.x functions (`plot_buffer`, `plot_buffer_traces` in `inc/plot.h`) still build as
`//:signal_sniper_plot` and will be removed later.

---

## Repository layout

```
inc/ssp/          Public C++ API: plot.h (plot, raster, Session, options), types.h (Signal, ...)
src/
  core/           Data access: typed Signal dispatch, sample scans, the min/max pyramid (lod)
  render/         Framebuffer, font, PNG, axes/ticks, colormaps, 1-D reduce and paint
  app/            Interactive app: content (xplot TraceContent, xraster RasterContent), screen
                  state, controller (input → actions), menu, overlays, cuts, render worker,
                  event loop (session)
  platform/       Window-system backend interface and the X11 implementation
  api/            plot / raster / save_png entry points
pybind_src/       Python module (signal_sniper_plot_py)
python/           Python tests (test_ssp.py) and demos (test_plot.py, ofdm_demo.py)
tests/ssp/        C++ unit, golden-image, event-loop and X11 tests
tests/golden/     Reference images for the golden tests
tools/            Benchmark, C++ demos, release build script, font/colormap generators
bazel/            Host X11 for Bazel; pinned pip requirements
docs/             This file, the PyPI readme, dependencies, architecture notes
inc/*.h, src/plot*.cc, src/trace_utils.cc, tests/test_plot.cc   1.x implementation (legacy)
```

[ARCHITECTURE_PROPOSAL.md](ARCHITECTURE_PROPOSAL.md) explains the design: zero-copy `Signal`
views, the reduce (worker) / paint (UI) split, the `Content` interface shared by xplot and
xraster, and the cut screen stack.

---

## Prerequisites

| Dependency | Purpose | Install (Ubuntu / WSL) |
|---|---|---|
| **Bazelisk** | Bazel version manager; reads `.bazeliskrc` (Bazel 9.0.0) | see below |
| **libx11-dev** | X11 headers and library | `sudo apt install libx11-dev` |
| **Docker** | Only for building release wheels | [docs.docker.com](https://docs.docker.com/engine/install/) |

Everything else (GoogleTest, pybind11, zlib, the Python 3.12 toolchain, numpy) is fetched
by Bazel; see [DEPENDENCIES.md](DEPENDENCIES.md).

```sh
curl -fsSL -o /usr/local/bin/bazel \
  https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64
chmod +x /usr/local/bin/bazel
```

---

## Building

```sh
bazel build //...                            # everything
bazel build //:ssp                           # C++ library (X11 window + headless)
bazel build //:ssp_core                      # headless core: no X11 dependency
bazel build //:signal_sniper_plot_py         # Python extension .so
bazel build -c opt //:signal_sniper_plot_wheel   # wheel for this machine only (see Releasing)
```

Use `-c opt` for anything performance-related; Bazel's default build is unoptimised and
4–5× slower.

---

## Testing

None of these need a display:

```sh
bazel test //:ssp_core_test //:ssp_golden_test //:ssp_py_test
```

- `ssp_core_test`: unit tests, plus the real event loop driven by a fake window (zoom,
  stop/resume, unzoom cache, menus, cuts).
- `ssp_golden_test`: renders known scenes and compares pixel for pixel with `tests/golden/`.
  After an intended visual change, inspect and accept the new images with
  `SSP_UPDATE_GOLDENS=1 bazel run //:ssp_golden_test` (on failure, the actual image is in
  `bazel-testlogs/ssp_golden_test/test.outputs/`).
- `ssp_py_test`: the Python module, headless.

The X11 backend test needs a display (`Xvfb :99 &` works):

```sh
DISPLAY=:99 bazel test //:ssp_x11_test --test_env=DISPLAY
```

The worker thread and UI share results through release/acquire counters; check changes there
with ThreadSanitizer (`setarch -R` works around TSan's trouble with high ASLR entropy):

```sh
bazel test //:ssp_core_test --copt=-fsanitize=thread --linkopt=-fsanitize=thread --copt=-O1 \
  --run_under="setarch x86_64 -R"
```

---

## Demos and benchmark

```sh
bazel run -c opt //:ssp_demo              # three 50M-sample traces (C++)
bazel run -c opt //:ssp_raster_demo       # 20000 x 4096 spectrogram (C++)
bazel run //:ofdm_demo                    # OFDM acquisition + tracking (Python); -- --png DIR to save
bazel run //:test_plot_py                 # short Python examples
bazel run -c opt //:ssp_bench             # render timings on 200M complex samples, no display
```

---

## Python dependencies

Declared in `bazel/python_deps/requirements.in`, pinned in
`bazel/python_deps/requirements_lock.txt`. After editing the `.in` file:

```sh
bazel run //:requirements.update
```

---

## Releasing

Release wheels are built inside PyPA's `manylinux_2_28` image, so they install on any x86_64
Linux with glibc 2.28 or newer (RHEL / CentOS Stream 8+, Ubuntu 20.04+). A wheel built
directly on a newer host would need that host's glibc.

```sh
tools/build_manylinux_wheel.sh   # needs Docker; writes dist/*.whl
```

Pushing a `v*` tag runs `.github/workflows/build.yml`, which runs the same script, checks the
wheel with `twine check`, and uploads it to TestPyPI. The `TEST_PYPI_API_TOKEN` secret must be
set in the repository's GitHub Actions settings. The version is set in `MODULE.bazel` and in
the `signal_sniper_plot_wheel` rule in `BUILD.bazel`.

---

## License

MIT — see [LICENSE](LICENSE).

## Contributing

Issues and pull requests are welcome. Please open an issue to discuss significant changes
before submitting a PR.

## Contact

[aj_horton@hotmail.com](mailto:aj_horton@hotmail.com)
