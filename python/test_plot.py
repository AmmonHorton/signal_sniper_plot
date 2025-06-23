import numpy as np
from signal_sniper_plot_py import plot_buffer

# Just an array, all other defaults
plot_buffer(np.array(list(range(100))), y_range=[0, 100])

# Complex input with optional args
plot_buffer(np.random.randn(1024).astype(np.complex128),
            xdelta=0.1, plot_title="IQ Sample", y_range=None, x_range=[0.5, 51.2])
