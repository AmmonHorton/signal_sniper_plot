# DSP Library Dependencies

This document lists the dependencies required to build and run the DSP Library.  
Most Python dependencies (**NumPy, SciPy, pytest**) are installed automatically  
when running the debug setup script (`scripts/setup_debug_env.sh`).  

Below are the dependencies that must be installed **manually** before building.

---

## ** Required System Dependencies**
These dependencies **are NOT handled by the virtual environment** and must be installed manually.

| Dependency  | Purpose |
|------------|---------|
| **GoogleTest** | Unit testing framework |
| **CMake** | Required for building the project |
| **ibx11-dev** | X11 Graphices library for rendering|
| **Pybind11** | C++ bindings for Python |


---

## ** Manual Installation Instructions**

### ** Ubuntu (WSL)**
Run the following command to install the required C++ libraries:

```sh
sudo apt update && sudo apt install -y \
    build-essential cmake python3-dev \
    libx11-dev libgtest-dev pybind11-dev
```

