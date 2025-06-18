cmake_minimum_required(VERSION 3.16)
project(Signal_Sniper_Plot VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options to enable/disable components
option(ENABLE_TESTS "Enable GoogleTest Unit Tests" OFF)
option(ENABLE_PYBIND "Enable Pybind11 Python Binding Compilation" OFF)
option(ENABLE_PYTHON_INSTALL "Enable Pybind11 python installation - Must have ENABLE_PYBIND set" OFF)
option(BUILD_RPM "Enable RPM packaging" OFF)
option(BUILD_DEB "Build a Debian package" OFF)
option(BUILD_WHEEL "Build a/install for python whele" OFF)


# Optimization and Warning Options
option(ENABLE_OPTIMIZED_BUILD "Enable -O2 optimization" OFF)
option(ENABLE_WARNINGS "Enable strict warnings (-Wall -Wextra)" OFF)

if (ENABLE_WARNINGS)
    message(STATUS "Enabling warnings for all compilation: -Wall -Wextra")
    add_compile_options(signal_sniper_plot PRIVATE -Wall -Wextra)
endif()

# Find X11 (for GUI)
find_package(X11 REQUIRED)

# Automatically find all source files
file(GLOB_RECURSE PLOT_SOURCES LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} CONFIGURE_DEPENDS src/*.cc)

# Debugging output: Show all detected sources
message(STATUS "Discovered source files: ${PLOT_SOURCES}")

# Core Signal Sniper Plot Library
add_library(signal_sniper_plot SHARED ${PLOT_SOURCES})

if (ENABLE_OPTIMIZED_BUILD)
    target_compile_options(signal_sniper_plot PRIVATE -O2)
endif()

target_include_directories(signal_sniper_plot PUBLIC inc)
target_link_libraries(signal_sniper_plot PRIVATE X11)

# GoogleTest (optional)
if(ENABLE_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    
    file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} CONFIGURE_DEPENDS tests/*.cc)
    message(STATUS "Discovered test source files: ${TEST_SOURCES}")

    if(TEST_SOURCES)
        add_executable(test_signal_sniper_plot ${TEST_SOURCES})
        target_link_libraries(test_signal_sniper_plot PRIVATE GTest::GTest GTest::Main signal_sniper_plot)
        add_test(NAME SignalSniperTests COMMAND test_signal_sniper_plot)
    else()
        message(WARNING "No test files found in tests/ directory. Skipping test target.")
    endif()
endif()

# Pybind11 Python Bindings (optional)
if(ENABLE_PYBIND)
    find_package(pybind11 REQUIRED)

    file(GLOB_RECURSE PYBIND_SOURCES LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} CONFIGURE_DEPENDS pybind_src/*.cc)
    message(STATUS "Discovered pybind source files: ${PYBIND_SOURCES}")

    # Add the pybind module and link against the other built library
    pybind11_add_module(signal_sniper_plot_py ${PYBIND_SOURCES})
    
    # Set RPATH so that signal_sniper_plot_py.so can find libsignal_sniper.so
    set_target_properties(signal_sniper_plot_py PROPERTIES INSTALL_RPATH "\$ORIGIN")
    target_link_libraries(signal_sniper_plot_py PRIVATE signal_sniper_plot)
    
    if (ENABLE_PYTHON_INSTALL)
        # Check if a virtual environment exists and use it if available
        if(EXISTS "${CMAKE_SOURCE_DIR}/.venv/bin/python3")
            set(Python3_EXECUTABLE "${CMAKE_SOURCE_DIR}/.venv/bin/python3")
            message(STATUS "Using virtual environment Python: ${Python3_EXECUTABLE}")
        endif()
        
        # Run find package which will set python3_SITELIB
        find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

        # Ensure Python3_SITELIB is set
        if(NOT Python3_SITELIB)
            message(FATAL_ERROR "CMake did not detect Python3_SITELIB. Ensure Python3 is properly found.")
        endif()
        message("Installing pybind modules at: ${Python3_SITELIB}")
        install(TARGETS signal_sniper_plot DESTINATION ${Python3_SITELIB})
        install(TARGETS signal_sniper_plot_py DESTINATION ${Python3_SITELIB})
    elseif (BUILD_WHEEL)
        # If building for python .whl, install to the current directory
        install(TARGETS signal_sniper_plot DESTINATION .)
        install(TARGETS signal_sniper_plot_py DESTINATION .)
    endif()
endif()

# RPM Packaging (optional)
if(BUILD_RPM)
    add_subdirectory(rpm)
endif()

# Debian Packaging (optional)
if(BUILD_DEB)
    add_subdirectory(deb)
endif()