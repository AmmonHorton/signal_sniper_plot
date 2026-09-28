# signal_sniper_plot_py

Fast interactive DSP plots for NumPy, in the spirit of XMidas `xplot` / `xraster` and SigPlot.
Pass an array, get a zoomable X11 window. Arrays of any size are read in place (no copies,
`np.memmap` works); the first view of a huge array renders progressively and every view after
that takes milliseconds.

> **Requires an X11 display** (`DISPLAY` must be set): Linux, WSL2 (WSLg or an X server such
> as VcXsrv), or `ssh -X`. `save_png` / `save_raster_png` need no display.

---

## Installation

```sh
pip install signal_sniper_plot_py
```

Linux x86_64 with glibc 2.28 or newer (RHEL / CentOS Stream / Alma / Rocky 8+, Ubuntu 20.04+,
Debian 10+), Python 3.12. `libX11` must be installed (it is on any desktop Linux;
`sudo apt install libx11-6` or `sudo dnf install libX11` otherwise).

---

## Quick start

```python
import numpy as np
import signal_sniper_plot_py as ssp

fs = 1e6
t = np.arange(1_000_000) / fs
iq = np.exp(2j * np.pi * 1e3 * t).astype(np.complex64)

ssp.plot(iq, xdelta=1 / fs, title="1 kHz tone")          # line plot (xplot)

spec = np.abs(np.fft.fft(iq.reshape(1000, 1000), axis=1))
ssp.raster(spec, cmode="20log", title="spectrogram")      # image (xraster)
```

Each call opens a window and returns when it is closed (Ctrl-C closes it too).

---

## Line plots: `plot`

```python
ssp.plot(*data, title="", cmode="auto", xstart=0.0, xdelta=1.0, names=None,
         xrange=None, yrange=None, index=False, thickness=1, grid=True, legend=True,
         phunits="rad", width=1000, height=600)
```

`data` is any number of:

- 1-D arrays — one trace each;
- 2-D arrays — one trace per row;
- `ssp.Signal` objects — a trace with its own x axis and style;
- lists of the above.

`xstart`/`xdelta` apply to bare arrays. A `Signal` carries its own:

```python
ssp.plot(ssp.Signal(rx, xstart=0.0, xdelta=1 / fs, name="rx"),
         ssp.Signal(tx, xstart=2e-3, xdelta=1 / fs, name="tx", style="dots", color="#ff8000"),
         title="rx vs tx", cmode="real", yrange=(-1.5, 1.5))
```

| Option | Meaning |
|---|---|
| `cmode` | `auto` (real data: `real`, complex: `mag`), `mag`, `phase`, `real`, `imag`, `ir` (imag vs real), `10log`, `20log` |
| `xrange`, `yrange` | Fixed `(lo, hi)` view; autoscaled when `None` |
| `index` | x axis in sample numbers instead of `xstart + i*xdelta` |
| `names` | Legend names for the traces, in order |
| `phunits` | Phase units: `rad`, `deg`, `cycles` |

`Signal(data, xstart=0.0, xdelta=1.0, name="", style="lines", color=None, thickness=0, visible=True)`
— `style` is `lines`, `dots` or `both`; `color` is `0xRRGGBB` or `"#RRGGBB"`.

---

## Rasters: `raster`

```python
ssp.raster(data, subsize=None, xstart=0.0, xdelta=1.0, ystart=0.0, ydelta=1.0, title="",
           cmode="auto", xrange=None, yrange=None, zrange=None, cmap="ramp", reduce="max",
           index=False, grid=False, phunits="rad", width=1000, height=600)
```

- `data` is a 2-D array with one frame per row, or a 1-D array with `subsize=` samples per
  frame. Strided views (`a[::2, 10:500]`) are read in place.
- Frame 0 is drawn at the top. Columns use `xstart`/`xdelta`; frames use `ystart`/`ydelta`.

| Option | Meaning |
|---|---|
| `zrange` | Fixed colour range `(lo, hi)`; autoscaled when `None` |
| `cmap` | `greyscale`, `ramp`, `colorwheel`, `spectrum`, `calewhite`, `hotdesat`, `sunset`, `hot`, `cold` |
| `reduce` | How a pixel combines the samples under it when zoomed out: `max`, `min`, `mean`, `maxabs`, `first` |

**Cuts:** in a raster window, `x` / `y` open a line plot of the row / column under the
pointer, limited to the current zoom box. It behaves like any `plot` window, and nothing done
in it changes the raster. `PgUp`/`PgDn` step to the neighbouring row/column; `Esc` returns.

---

## Without a display: `save_png`, `save_raster_png`

```python
ssp.save_png("tone.png", iq, xdelta=1 / fs, title="1 kHz tone")
ssp.save_raster_png("spec.png", spec, cmode="20log")
```

Same arguments as `plot` / `raster`, after the output path.

---

## In the window

Press `?` for the full list, or middle-click (or `m`) for the menu.

| Input | Action |
|---|---|
| Left drag | Zoom to the box (right-click steps back out; `Home` unzooms fully) |
| Wheel | Zoom x around the pointer (Shift: y, Ctrl: both) |
| Shift + left drag, arrows | Pan |
| Left click | Set a marker; the readout then shows dx / dy from it |
| `1` … `7` | Magnitude, Phase, Real, Imag, Imag vs Real, 10·log10, 20·log10 |
| Legend click | Show / hide a trace (right-click: lines → dots → both) |
| `l`, `g`, `i`, `a` | Legend, grid, index x axis, readout x / index / 1/x |
| `c`, `[` `]` | Raster colormap; slide the colour range |
| Any click, `Esc` / `Space` | Stop / resume a render in progress |
| `Ctrl-S`, `q` | Save PNG, quit |

---

## Data types

Any of `int8`, `uint8`, `int16`, `uint16`, `int32`, `uint32`, `int64`, `uint64`, `float32`,
`float64`, `complex64`, `complex128` and `bool`, in any byte order.
`float16` and `complex256` are rejected.

---

## Earlier API

`plot_buffer(data, xstart, xdelta, plot_title, line_thickness, y_range, x_range, num_traces)`
from 1.x still works and opens the new window; `plot` is the replacement.

---

## License

MIT
