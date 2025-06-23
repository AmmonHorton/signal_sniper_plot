// trace_utils.cpp

#include "trace_utils.h"
#include <stdexcept>
#include <cstdint>
#include <cmath>

namespace xplot {

Trace process_buffer(const void* data,
                     std::size_t num_elements,
                     std::size_t element_size_bytes,
                     bool is_complex,
                     bool is_float,
                     double xstart,
                     double xdelta) {
    Trace trace;
    trace.xstart = xstart;
    trace.xdelta = xdelta;
    trace.samples.resize(num_elements);

    if (is_float) {
        if (element_size_bytes == 4 || (is_complex && element_size_bytes == 8)) {
            const float* fdata = reinterpret_cast<const float*>(data);
            for (std::size_t i = 0; i < num_elements; ++i) {
                float re = fdata[is_complex ? 2 * i : i];
                float im = is_complex ? fdata[2 * i + 1] : 0.0f;
                trace.samples[i] = {xstart + i * xdelta, re, im};
            }
        } else if (element_size_bytes == 8 || (is_complex && element_size_bytes == 16)) {
            const double* ddata = reinterpret_cast<const double*>(data);
            for (std::size_t i = 0; i < num_elements; ++i) {
                float re = static_cast<float>(ddata[is_complex ? 2 * i : i]);
                float im = is_complex ? static_cast<float>(ddata[2 * i + 1]) : 0.0f;
                trace.samples[i] = {xstart + i * xdelta, re, im};
            }
        } else {
            throw std::runtime_error("Unsupported float element size: " + std::to_string(element_size_bytes));
        }
    } else {
        if (element_size_bytes == 2 || (is_complex && element_size_bytes == 4)) {
            const int16_t* sdata = reinterpret_cast<const int16_t*>(data);
            for (std::size_t i = 0; i < num_elements; ++i) {
                float re = static_cast<float>(sdata[is_complex ? 2 * i : i]);
                float im = is_complex ? static_cast<float>(sdata[2 * i + 1]) : 0.0f;
                trace.samples[i] = {xstart + i * xdelta, re, im};
            }
        } else if (element_size_bytes == 4 || (is_complex && element_size_bytes == 8)) {
            const int32_t* idata = reinterpret_cast<const int32_t*>(data);
            for (std::size_t i = 0; i < num_elements; ++i) {
                float re = static_cast<float>(idata[is_complex ? 2 * i : i]);
                float im = is_complex ? static_cast<float>(idata[2 * i + 1]) : 0.0f;
                trace.samples[i] = {xstart + i * xdelta, re, im};
            }
        } else if (element_size_bytes == 8 || (is_complex && element_size_bytes == 16)) {
            const int64_t* idata = reinterpret_cast<const int64_t*>(data);
            for (std::size_t i = 0; i < num_elements; ++i) {
                float re = static_cast<float>(idata[is_complex ? 2 * i : i]);
                float im = is_complex ? static_cast<float>(idata[2 * i + 1]) : 0.0f;
                trace.samples[i] = {xstart + i * xdelta, re, im};
            }
        } else {
            throw std::runtime_error("Unsupported int element size: " + std::to_string(element_size_bytes));
        }
    }

    return trace;
}

}  // namespace xplot
