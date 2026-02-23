# Signal Sniper Plot

An interactive X11 DSP plotting library for C++ and Python, inspired by XMidas `xplot`. Designed for real-time signal debugging: open a window, stream traces left-to-right as data arrives, zoom/pan freely, and toggle trace visibility from the legend — all without blocking your application.

---

## Features

- **Streaming render** — traces draw column-by-column so the window is responsive immediately; a per-trace decimation cache makes subsequent renders instant
- **Zoom / pan** — left-drag to zoom in (up to 5 levels), right-click to step back out
- **Plot modes** — Magnitude, Real, Imaginary, Phase, and I/Q scatter, switchable at runtime from the toolbar
- **Legend** — click a swatch to hide/show a trace; right-click to toggle LINES ↔ DOTS style
- **Crosshair readout** — cursor coordinates shown live in the toolbar
- **Python bindings** — full pybind11 wrapper; accepts NumPy arrays directly

---

## Repository layout

```
inc/            Public C++ headers (plot.h, plot_session.h, trace_utils.h)
src/            Library implementation (plot.cc, plot_session.cc, trace_utils.cc)
pybind_src/     Pybind11 Python extension (module.cc)
tests/          GoogleTest unit/integration tests
python/         Python demo script (test_plot.py)
bazel/
  system_libs/  BUILD file exposing the host X11 installation to Bazel
  python_deps/  pip requirements files
docs/           This file and DEPENDENCIES.md
.github/
  workflows/    CI: builds and publishes the PyPI wheel on version tags
```

---

## Prerequisites

| Dependency | Purpose | Install (Ubuntu/WSL) |
|---|---|---|
| **Bazelisk** | Bazel version manager; reads `.bazeliskrc` | see below |
| **libx11-dev** | X11 headers and shared library | `sudo apt install libx11-dev` |

All other dependencies (GoogleTest, pybind11, rules_python, etc.) are fetched automatically by Bazel from the Bazel Central Registry.

### Install Bazelisk

```sh
curl -fsSL -o /usr/local/bin/bazel \
  https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64
chmod +x /usr/local/bin/bazel
```

Bazelisk reads `.bazeliskrc` at the repo root and downloads the pinned Bazel version (9.0.0) on first use.

---

## Building

```sh
# Build the full C++ library and Python extension
bazel build //...

# Build only the C++ library
bazel build //:signal_sniper_plot

# Build only the Python extension .so
bazel build //:signal_sniper_plot_py

# Build the distributable Python wheel
bazel build //:signal_sniper_plot_wheel
```

---

## Running tests

The test opens an interactive X11 window — a live `DISPLAY` is required.

```sh
bazel test //:test_signal_sniper_plot --test_output=streamed
```

---

## Running the Python demo

```sh
bazel run //:test_plot_py
```

---

## Managing Python dependencies

pip dependencies are declared in `bazel/python_deps/requirements.in` and pinned in `bazel/python_deps/requirements_lock.txt`. To regenerate the lock file after editing the `.in` file:

```sh
bazel run //:requirements.update
```

---

## Publishing a release

Pushing a tag of the form `v*` triggers the GitHub Actions workflow (`.github/workflows/build.yml`), which:

1. Installs system deps and Bazelisk on the runner
2. Builds the wheel with `bazel build //:signal_sniper_plot_wheel`
3. Runs `auditwheel repair` to apply the correct `manylinux` platform tag
4. Checks the wheel with `twine check`
5. Uploads to TestPyPI via `twine upload`

The `TEST_PYPI_API_TOKEN` secret must be set in the repository's GitHub Actions settings.

---

## License

MIT — see [LICENSE](LICENSE).

## Contributing

Issues and pull requests are welcome. Please open an issue to discuss significant changes before submitting a PR.

## Contact

[aj_horton@hotmail.com](mailto:aj_horton@hotmail.com)
