"""Deprecated alias: the module is now `signal_sniper_plot`. This alias goes away in 3.0."""

import warnings

from signal_sniper_plot import *  # noqa: F401,F403

warnings.warn("import signal_sniper_plot instead of signal_sniper_plot_py", DeprecationWarning, stacklevel=2)
