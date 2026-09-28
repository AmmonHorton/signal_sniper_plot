"""Manual demo (needs a display):  bazel run //:test_plot_py"""
import numpy as np

import signal_sniper_plot_py as ssp

fs = 1e6
t = np.arange(5_000_000) / fs
tone = np.exp(2j * np.pi * 1e3 * t).astype(np.complex64)
noise = (np.random.randn(t.size) * 0.1).astype(np.float32)

# One array, all defaults.
ssp.plot(np.arange(100), yrange=(0, 100))

# Several traces with their own axes and styles.
ssp.plot(ssp.Signal(tone, xdelta=1 / fs, name="tone"),
         ssp.Signal(noise, xstart=1.0, xdelta=1 / fs, name="noise", style="dots"),
         title="Two traces", cmode="real")

# The original API still works.
ssp.plot_buffer(np.random.randn(1024).astype(np.complex128),
                xdelta=0.1, plot_title="IQ Sample", x_range=(0.5, 51.2))
