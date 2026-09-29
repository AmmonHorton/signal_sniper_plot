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

Python (`pip install signal_sniper_plot_py`; see [python/README.md](python/README.md) for the full reference):

```python
import numpy as np
import signal_sniper_plot as ssp

ssp.plot(iq, xdelta=1 / fs, cmode="mag")
ssp.plot(ssp.Signal(rx, xdelta=1 / fs, name="rx"), ssp.Signal(tx, xstart=2e-3, xdelta=1 / fs, name="tx"))
ssp.raster(spectrogram, cmode="20log", cmap="ramp")
```

C++ (`#include "ssp/plot.h"`; Bazel target `//ssp`, or the RPM / .deb with pkg-config):

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
The 1.x C++ functions (`plot_buffer`, `plot_buffer_traces`) were removed in 2.0; Python's
`plot_buffer` still works.

---

## Repository layout

```
ssp/                  The C++ library. Public API: plot.h, types.h        //ssp, //ssp:headless
  core/               Typed sample access, scans, the min/max pyramid      //ssp/core
  render/             Framebuffer, font, PNG, axes, colormaps, 1-D paint    //ssp/render
  app/                xplot / xraster content, screens, input, menu,
                      overlays, cuts, render worker, event loop            //ssp/app
    testdata/golden/  Reference images for the golden test
  platform/           Window-system interface; x11/ implements it          //ssp/platform/x11
python/               Python module signal_sniper_plot, its test, PyPI readme, pip requirements
examples/             cpp/ (plot and raster demos), python/ (basic, OFDM)
benchmarks/           Rendering benchmark
packaging/            Wheel, RPM and .deb rules; build.sh builds them in containers
third_party/          Host libraries for Bazel (X11)
tools/                Generators for the embedded font and colormap tables
docs/                 Architecture and dependency notes
```

Tests live next to the code they test (`ssp/core/lod.cc` / `lod_test.cc`). Only `//ssp`,
`//ssp:headless` and the Python module are public; the rest is visible inside `ssp/` only.
[docs/architecture.md](docs/architecture.md) explains the design: zero-copy `Signal` views, the
reduce (worker) / paint (UI) split, the `Content` interface shared by xplot and xraster, and the
cut screen stack.

---

## Prerequisites

| Dependency | Purpose | Install (Ubuntu / WSL) |
|---|---|---|
| **Bazelisk** | Bazel version manager; reads `.bazeliskrc` (Bazel 9.0.0) | see below |
| **libx11-dev** | X11 headers and library | `sudo apt install libx11-dev` |
| **Docker** | Only for release builds (`packaging/build.sh`) | [docs.docker.com](https://docs.docker.com/engine/install/) |

Everything else (GoogleTest, pybind11, zlib, the Python 3.12 toolchain, numpy) is fetched
by Bazel; see [docs/dependencies.md](docs/dependencies.md).

```sh
curl -fsSL -o /usr/local/bin/bazel \
  https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64
chmod +x /usr/local/bin/bazel
```

---

## Building

```sh
bazel build //...                          # everything (except the RPM, which needs rpmbuild)
bazel build //ssp                          # C++ library: X11 windows + headless PNG
bazel build //ssp:headless                 # save_png / save_raster_png only, no X11
bazel build //python:signal_sniper_plot    # Python extension .so
```

C++17 and the other shared flags live in `.bazelrc`. Use `-c opt` for anything
performance-related; Bazel's default build is unoptimised and 4–5× slower.

---

## Testing

```sh
bazel test //...        # every test; none needs a display
```

- `//ssp/core`, `//ssp/render`: unit tests.
- `//ssp/app`: the controller, xraster and cuts, and the real event loop driven by a fake
  window (zoom, stop/resume, unzoom cache, menus, cuts).
- `//ssp/app:golden_test`: renders known scenes and compares them pixel for pixel with
  `ssp/app/testdata/golden/`. After an intended visual change, inspect and accept the new
  images with `SSP_UPDATE_GOLDENS=1 bazel run //ssp/app:golden_test` (on failure, the actual
  image is in `bazel-testlogs/ssp/app/golden_test/test.outputs/`).
- `//python:signal_sniper_plot_test`: the Python module, headless.

Sanitizers are configured in `.bazelrc`. The worker thread and UI share results through
release/acquire counters, so run ThreadSanitizer after touching that:

```sh
bazel test --config=tsan //ssp/...
bazel test --config=asan //ssp/...    # AddressSanitizer + UndefinedBehaviorSanitizer
```

The X11 backend test needs a display (`Xvfb :99 &` works):

```sh
DISPLAY=:99 bazel test //ssp/platform/x11:x11_backend_test --test_env=DISPLAY
```

---

## Examples and benchmark

```sh
bazel run -c opt //examples/cpp:plot_demo       # three 50M-sample traces
bazel run -c opt //examples/cpp:raster_demo     # 20000 x 4096 spectrogram
bazel run //examples/python:ofdm                # OFDM acquisition + tracking; -- --png DIR saves figures
bazel run //examples/python:basic               # short Python examples
bazel run -c opt //benchmarks:render_bench      # render timings on 200M complex samples, no display
```

---

## Python dependencies

Declared in `python/requirements.in`, pinned in `python/requirements_lock.txt`. After editing
the `.in` file:

```sh
bazel run //python:requirements.update
```

---

## Releasing

Everything distributable is defined in `packaging/`: the Python wheel, and RPM / .deb packages
of the C++ library (shared library, headers and a pkg-config file, so a program builds with
`g++ -std=c++17 app.cc $(pkg-config --cflags --libs signal_sniper_plot)`).

Each is built inside the oldest distribution it should install on, since a package built on a
newer host needs that host's glibc:

```sh
packaging/build.sh wheel   # manylinux_2_28: pip on any x86_64 Linux with glibc >= 2.28
packaging/build.sh rpm     # Rocky Linux 8: RHEL / CentOS Stream / Alma / Rocky 8 and newer
packaging/build.sh deb     # Ubuntu 20.04 and newer
packaging/build.sh all     # needs Docker; writes dist/
```

Pushing a `v*` tag runs `.github/workflows/build.yml`: it runs `packaging/build.sh all`, checks
the wheel with `twine check` and uploads it to TestPyPI, and attaches the RPM and .deb to the
workflow run. The `TEST_PYPI_API_TOKEN` secret must be set in the repository's GitHub Actions
settings. The version is set in `MODULE.bazel` and `packaging/BUILD.bazel`.

---

## License

MIT — see [LICENSE](LICENSE).

## Contributing

Issues and pull requests are welcome. Please open an issue to discuss significant changes
before submitting a PR.

## Contact

[aj_horton@hotmail.com](mailto:aj_horton@hotmail.com)
