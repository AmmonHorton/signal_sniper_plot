# signal_sniper_plot_py

An interactive X11 DSP plotting library for Python. Pass a NumPy array, get a live, zoomable signal window — designed for real-time debugging of DSP pipelines.

Backed by a C++ rendering engine with per-trace streaming, a decimation cache, and full zoom/pan/mode switching from the window UI.

> **Requires a live X11 display** (`DISPLAY` must be set). Works natively on Linux and under WSL2 with an X server (e.g. VcXsrv, X410).

---

## Installation

```sh
pip install signal_sniper_plot_py
```

**System requirement:** `libX11` must be present on the host. On Ubuntu/Debian:

```sh
sudo apt install libx11-dev
```

---

## Quick start

```python
import numpy as np
from signal_sniper_plot_py import plot_buffer

# Plot a simple sine wave
t = np.linspace(0, 1, 4096, dtype=np.float32)
plot_buffer(np.sin(2 * np.pi * 50 * t), xdelta=1/4096, plot_title="50 Hz sine")
```

```python
# Plot a complex IQ signal
iq = (np.random.randn(1024) + 1j * np.random.randn(1024)).astype(np.complex64)
plot_buffer(iq, xdelta=0.001, plot_title="IQ noise")
```

```python
# Two interleaved traces in one array
data = np.concatenate([signal_a, signal_b])   # shape: (2*N,)
plot_buffer(data, num_traces=2, plot_title="A vs B")
```

---

## API reference

### `plot_buffer`

```python
plot_buffer(
    data,
    xstart=0.0,
    xdelta=1.0,
    plot_title="Plot",
    line_thickness=2,
    y_range=None,
    x_range=None,
    num_traces=1,
)
```

Opens a blocking interactive plot window. Returns when the window is closed.

#### Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | `np.ndarray` | — | 1-D C-contiguous array of samples. See supported dtypes below. |
| `xstart` | `float` | `0.0` | X-axis value assigned to the first sample. |
| `xdelta` | `float` | `1.0` | X-axis increment between consecutive samples (e.g. `1/sample_rate`). |
| `plot_title` | `str` | `"Plot"` | Text shown in the window title bar and above the plot area. |
| `line_thickness` | `int` | `2` | Pixel thickness for line traces; dot radius for DOTS-style traces. |
| `y_range` | `(float, float) \| None` | `None` | Fix the Y-axis to `(min, max)`. Auto-scales from data when `None`. |
| `x_range` | `(float, float) \| None` | `None` | Fix the X-axis to `(min, max)`. Auto-scales from data when `None`. |
| `num_traces` | `int` | `1` | Number of traces interleaved in `data`. The array is split evenly: trace 0 occupies `data[0:N]`, trace 1 `data[N:2N]`, etc. |

#### Supported dtypes

| NumPy dtype | Interpreted as |
|---|---|
| `float32` | Real samples (32-bit float) |
| `float64` | Real samples (64-bit double) |
| `int16` | Real samples (16-bit integer) |
| `int32` | Real samples (32-bit integer) |
| `int64` | Real samples (64-bit integer) |
| `complex64` | Interleaved real + imag (2 × float32) |
| `complex128` | Interleaved real + imag (2 × float64) |

---

## Window controls

Once the plot window is open:

| Action | Effect |
|---|---|
| Left-drag | Zoom into the selected rectangle (up to 5 levels) |
| Right-click (plot area) | Step back out one zoom level |
| Middle-click | Pause / resume the streaming render |
| Left-click legend swatch | Toggle trace visibility |
| Right-click legend swatch | Toggle trace style: LINES ↔ DOTS |

### Toolbar buttons

| Button | Effect |
|---|---|
| **Save PNG** | Saves the current view to `plot_out.ppm` |
| **Cycle X-Axis** | Toggle x-axis labels between time values and sample indices |
| **Magnitude** | Display √(re² + im²) |
| **Real** | Display real component only |
| **Imag** | Display imaginary component only |
| **Phase** | Display atan2(re, im) in radians |
| **Imag vs Real** | I/Q scatter plot |

---

## Examples

### Fixed axes with a complex signal

```python
import numpy as np
from signal_sniper_plot_py import plot_buffer

fs = 10_000          # sample rate, Hz
N  = 50_000
t  = np.arange(N) / fs

# Chirp: frequency sweeps from 100 Hz to 2 kHz
f  = np.linspace(100, 2000, N)
iq = np.exp(1j * 2 * np.pi * np.cumsum(f) / fs).astype(np.complex64)

plot_buffer(
    iq,
    xstart=0.0,
    xdelta=1/fs,
    plot_title="Chirp",
    y_range=(-1.1, 1.1),
)
```

### Comparing two signals

```python
import numpy as np
from signal_sniper_plot_py import plot_buffer

N     = 8192
t     = np.linspace(0, 1, N, dtype=np.float32)
clean = np.sin(2 * np.pi * 440 * t)
noisy = clean + 0.3 * np.random.randn(N).astype(np.float32)

# Interleave: [clean_0, clean_1, ..., noisy_0, noisy_1, ...]
plot_buffer(
    np.concatenate([clean, noisy]),
    num_traces=2,
    xdelta=1/N,
    plot_title="Clean vs Noisy",
    y_range=(-2.0, 2.0),
)
```

---

## Requirements

- Python >= 3.10
- NumPy
- Linux with `libX11` and a running X11 display (`$DISPLAY`)

---

## License

MIT
