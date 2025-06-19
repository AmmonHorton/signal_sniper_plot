import numpy as np
from signal_sniper_plot_py import plot_buffer

# Just an array, all other defaults
plot_buffer(np.random.randn(2048).astype(np.int16))

# Complex input with optional args
plot_buffer(np.random.randn(1024).astype(np.float32),
            xdelta=0.001, plot_title="IQ Sample", y_range=None)
