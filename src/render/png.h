/// @file png.h
/// @brief Minimal PNG I/O for Framebuffer (8-bit RGB).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "render/framebuffer.h"

namespace ssp {

std::vector<uint8_t> encode_png(const Framebuffer& fb);

/// @throws std::runtime_error if the file can't be written.
void write_png(const Framebuffer& fb, const std::string& path);

/// @brief Decode PNGs as written by encode_png (8-bit RGB, no interlace, filter 0).
/// @throws std::runtime_error on anything else.
Framebuffer decode_png(const std::vector<uint8_t>& bytes);
Framebuffer read_png(const std::string& path);

}  // namespace ssp
