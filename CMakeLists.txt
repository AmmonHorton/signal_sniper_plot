cmake_minimum_required(VERSION 3.12)
project(xplotter LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories(inc)

add_library(xplot STATIC src/plot.cpp)

add_executable(test_plot test/main.cpp)
target_link_libraries(test_plot xplot X11)