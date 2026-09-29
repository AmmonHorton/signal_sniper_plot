"""Headless tests for the Python bindings (no display needed)."""
import os
import tempfile
import unittest

import numpy as np

import signal_sniper_plot as ssp


class SavePng(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.n = 0

    def png(self, *data, **kw):
        """Render and return the PNG bytes."""
        self.n += 1
        path = os.path.join(self.dir, f"{self.n}.png")
        ssp.save_png(path, *data, width=400, height=300, **kw)
        with open(path, "rb") as f:
            return f.read()

    def test_every_supported_dtype(self):
        base = (np.arange(1000) % 100).astype(np.float64)
        for dt in ["i1", "u1", "i2", "u2", "i4", "u4", "i8", "u8", "f4", "f8", "c8", "c16", "?"]:
            with self.subTest(dtype=dt):
                self.assertTrue(self.png(base.astype(dt)).startswith(b"\x89PNG"))

    def test_unsigned_values_are_not_sign_flipped(self):
        # Regression: uint16 used to be read as int16.
        v = np.linspace(0, 65535, 5000)
        u = v.astype(np.uint16)
        self.assertEqual(self.png(u), self.png(u.astype(np.float64)))

    def test_strided_and_reversed_views_match_copies(self):
        v = np.sin(np.arange(30000) * 0.01)
        self.assertEqual(self.png(v[::3]), self.png(np.ascontiguousarray(v[::3])))
        self.assertEqual(self.png(v[::-1]), self.png(v[::-1].copy()))

    def test_big_endian_matches_native(self):
        v = np.cos(np.arange(5000) * 0.02).astype(np.float32)
        self.assertEqual(self.png(v.astype(">f4")), self.png(v))

    def test_2d_is_one_trace_per_row_and_equals_separate_arrays(self):
        m = np.vstack([np.sin(np.arange(2000) * 0.01), np.cos(np.arange(2000) * 0.01)])
        self.assertEqual(self.png(m), self.png(m[0], m[1]))
        self.assertEqual(self.png(m[:, ::2]), self.png(m[0, ::2], m[1, ::2]))  # strided rows

    def test_signals_and_options(self):
        a = np.exp(1j * np.arange(4000) * 0.01).astype(np.complex64)
        b = ssp.Signal(np.arange(3000, dtype=np.int16), xstart=10, xdelta=0.5, name="ramp",
                       style="dots", color="#ff8000", thickness=2)
        out = self.png([a, b], title="t", cmode="20log", names=["tone", "ramp"],
                       xrange=(0, 100), yrange=(-80, 80), phunits="deg")
        self.assertTrue(out.startswith(b"\x89PNG"))

    def test_ir_mode(self):
        iq = np.exp(1j * (np.pi / 4 + np.pi / 2 * np.random.randint(0, 4, 50000))).astype(np.complex64)
        self.assertTrue(self.png(iq, cmode="ir", xrange=(1000, 20000)).startswith(b"\x89PNG"))

    def test_bad_input_raises_python_errors(self):
        v = np.zeros(10)
        with self.assertRaises(TypeError):
            self.png(v.astype(np.float16))
        with self.assertRaises(ValueError):
            self.png(np.zeros((2, 2, 2)))
        with self.assertRaises(ValueError):
            self.png()
        with self.assertRaises(ValueError):
            self.png(v, cmode="bogus")
        with self.assertRaises(ValueError):
            self.png(v, names=["a", "b"])
        with self.assertRaises(ValueError):
            self.png(v, yrange=(1, 1))
        with self.assertRaises(ValueError):
            ssp.Signal(v, style="zigzag")
        with self.assertRaises(ValueError):
            self.png(v, xdelta=0)


class Raster(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.n = 0

    def png(self, data, **kw):
        self.n += 1
        path = os.path.join(self.dir, f"r{self.n}.png")
        ssp.save_raster_png(path, data, width=400, height=300, **kw)
        with open(path, "rb") as f:
            return f.read()

    def test_2d_and_1d_with_subsize_match(self):
        m = np.outer(np.sin(np.arange(120) * 0.1), np.cos(np.arange(80) * 0.05)).astype(np.float32)
        self.assertEqual(self.png(m), self.png(m.ravel(), subsize=80))

    def test_views_read_in_place_match_copies(self):
        big = np.random.default_rng(0).standard_normal((200, 300))
        view = big[10:150:2, 5:250:3]  # strided rows and columns
        self.assertEqual(self.png(view), self.png(np.ascontiguousarray(view)))
        self.assertEqual(self.png(big[::-1]), self.png(big[::-1].copy()))

    def test_options(self):
        spec = (np.random.default_rng(1).standard_normal((64, 128)) +
                1j * np.random.default_rng(2).standard_normal((64, 128))).astype(np.complex64)
        out = self.png(spec, cmode="20log", cmap="hot", reduce="mean", zrange=(-40, 20),
                       xstart=-64, xdelta=1.0, ystart=0, ydelta=0.5, title="spec", grid=True)
        self.assertTrue(out.startswith(b"\x89PNG"))

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            self.png(np.zeros(100))  # 1-D without subsize
        with self.assertRaises(ValueError):
            self.png(np.zeros((4, 5)), subsize=3)
        with self.assertRaises(ValueError):
            self.png(np.zeros((4, 5)), cmap="rainbow")
        with self.assertRaises(ValueError):
            self.png(np.zeros((4, 5)), cmode="ir")


class Window(unittest.TestCase):
    def test_no_display_is_a_runtime_error_not_a_crash(self):
        saved = os.environ.pop("DISPLAY", None)
        try:
            with self.assertRaises(RuntimeError):
                ssp.plot(np.zeros(10))
            with self.assertRaises(RuntimeError):
                ssp.raster(np.zeros((4, 4)))
        finally:
            if saved is not None:
                os.environ["DISPLAY"] = saved

    def test_plot_buffer_shim_validates_before_opening(self):
        with self.assertRaises(ValueError):
            ssp.plot_buffer(np.zeros(10), num_traces=3)


class OldImportName(unittest.TestCase):
    def test_signal_sniper_plot_py_still_works_but_warns(self):
        import importlib
        with self.assertWarns(DeprecationWarning):
            import signal_sniper_plot_py
            importlib.reload(signal_sniper_plot_py)
        self.assertIs(signal_sniper_plot_py.plot, ssp.plot)
        self.assertIs(signal_sniper_plot_py.Signal, ssp.Signal)


if __name__ == "__main__":
    unittest.main()
