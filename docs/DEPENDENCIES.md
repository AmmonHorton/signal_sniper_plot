# Dependencies

All build-time dependencies except the two listed below are fetched automatically by Bazel from the Bazel Central Registry on first build. No manual installation of GoogleTest, pybind11, zlib, or Python toolchains is required.

At run time the library (and the Python wheel) needs only `libX11` and the C/C++ runtime; everything else is linked statically.

---

## Required system dependencies

These must be installed on the host before building.

| Dependency | Purpose | Install (Ubuntu / WSL) |
|---|---|---|
| **Bazelisk** | Bazel version manager; reads `.bazeliskrc` and downloads the correct Bazel release | See below |
| **libx11-dev** | X11 headers and shared library used by the rendering engine | `sudo apt install libx11-dev` |

### Installing Bazelisk

```sh
curl -fsSL -o /usr/local/bin/bazel \
  https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64
chmod +x /usr/local/bin/bazel
```

---

## Bazel-managed dependencies

These are declared in `MODULE.bazel` and require no manual action.

| Dependency | Version | Purpose |
|---|---|---|
| `rules_cc` | 0.2.16 | C++ build rules |
| `zlib` | 1.3.1 | PNG encoding (`save_png`, Ctrl-S). Built from source and linked statically: not needed where the library or wheel is used |
| `googletest` | 1.17.0 | Unit testing framework |
| `rules_python` | 1.8.3 | Python 3.12 hermetic toolchain and pip integration |
| `pybind11_bazel` | 3.0.0 | `pybind_extension()` rule for building the Python `.so` |

---

## Python pip dependencies

Declared in `bazel/python_deps/requirements.in` and pinned in `bazel/python_deps/requirements_lock.txt`. Managed entirely by Bazel; no virtualenv is needed.

| Package | Purpose |
|---|---|
| `numpy` | Required by the Python demo script and the pybind11 numpy interface |

To regenerate the lock file after editing `requirements.in`:

```sh
bazel run //:requirements.update
```
