/// @file dispatch.h
/// @brief The one place that turns a run-time DType into a compile-time element type.
#pragma once

#include <cstdint>
#include <stdexcept>

#include "ssp/types.h"

namespace ssp {

template <class T> struct Tag { using type = T; };

/// @brief Call `f(Tag<T>{})` with T matching `dt`.
template <class F>
decltype(auto) dispatch(DType dt, F&& f) {
    switch (dt) {
        case DType::I8:  return f(Tag<int8_t>{});
        case DType::U8:  return f(Tag<uint8_t>{});
        case DType::I16: return f(Tag<int16_t>{});
        case DType::U16: return f(Tag<uint16_t>{});
        case DType::I32: return f(Tag<int32_t>{});
        case DType::U32: return f(Tag<uint32_t>{});
        case DType::I64: return f(Tag<int64_t>{});
        case DType::U64: return f(Tag<uint64_t>{});
        case DType::F32: return f(Tag<float>{});
        case DType::F64: return f(Tag<double>{});
    }
    throw std::logic_error("unknown DType");
}

/// @brief Bytes per element (per real/imag half for complex data).
inline std::size_t dtype_size(DType dt) {
    return dispatch(dt, [](auto tag) { return sizeof(typename decltype(tag)::type); });
}

}  // namespace ssp
