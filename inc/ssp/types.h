/// @file types.h
/// @brief Public value types shared by xplot and xraster: Signal, CMode, Style, Range.
#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ssp {

/// @brief Element type of a sample (or of each real/imag half of a complex sample).
enum class DType : uint8_t { I8, U8, I16, U16, I32, U32, I64, U64, F32, F64 };

/// @brief Which view of the data is plotted. Two-letter names follow XMidas/SigPlot.
enum class CMode : uint8_t {
    Auto,   ///< Real for real data, Mag if any signal is complex.
    Mag,    ///< "Ma" |x|
    Phase,  ///< "Ph" arg(x)
    Real,   ///< "Re"
    Imag,   ///< "Im"
    IR,     ///< "IR" imaginary vs real (not implemented yet)
    Log10,  ///< "Lo" 10*log10|x|
    Log20,  ///< "L2" 20*log10|x|
};

enum class PhaseUnits : uint8_t { Radians, Degrees, Cycles };

/// @brief How a trace is drawn.
enum class Style : uint8_t { Lines, Dots, LinesDots };

/// @brief A closed interval [lo, hi] in data units.
struct Range {
    double lo = 0.0;
    double hi = 0.0;
};

template <class T> struct dtype_of;
template <> struct dtype_of<int8_t>   { static constexpr DType value = DType::I8; };
template <> struct dtype_of<uint8_t>  { static constexpr DType value = DType::U8; };
template <> struct dtype_of<int16_t>  { static constexpr DType value = DType::I16; };
template <> struct dtype_of<uint16_t> { static constexpr DType value = DType::U16; };
template <> struct dtype_of<int32_t>  { static constexpr DType value = DType::I32; };
template <> struct dtype_of<uint32_t> { static constexpr DType value = DType::U32; };
template <> struct dtype_of<int64_t>  { static constexpr DType value = DType::I64; };
template <> struct dtype_of<uint64_t> { static constexpr DType value = DType::U64; };
template <> struct dtype_of<float>    { static constexpr DType value = DType::F32; };
template <> struct dtype_of<double>   { static constexpr DType value = DType::F64; };

/// @brief A non-owning, typed view of evenly spaced samples plus how to draw them.
///
/// Samples are read in place: nothing is copied or converted. Sample i lives at element
/// `i * stride` (times 2 for complex, interleaved re/im) and has x = xstart + i * xdelta.
/// The caller keeps the memory alive for the lifetime of the plot, or hands ownership
/// to `keepalive`.
struct Signal {
    const void* data = nullptr;
    DType dtype = DType::F32;
    bool complex = false;
    std::size_t n = 0;          ///< Number of samples (complex pairs count once).
    std::ptrdiff_t stride = 1;  ///< Distance between samples, in samples (> 0).
    double xstart = 0.0;
    double xdelta = 1.0;        ///< Must be > 0.

    std::string name;           ///< Legend label; "Trace <i>" when empty.
    Style style = Style::Lines;
    bool visible = true;
    std::optional<uint32_t> color;  ///< 0xRRGGBB; palette colour when empty.
    int thickness = 0;              ///< Pixels; 0 = use the plot-wide thickness.

    std::shared_ptr<const void> keepalive;  ///< Optional owner of `data`.

    Signal() = default;

    /// @brief View of raw memory with a run-time element type.
    Signal(const void* p, DType dt, bool is_complex, std::size_t count,
           double x0 = 0.0, double dx = 1.0, std::string label = {})
        : data(p), dtype(dt), complex(is_complex), n(count), xstart(x0), xdelta(dx),
          name(std::move(label)) {}

    template <class T>
    Signal(const T* p, std::size_t count, double x0 = 0.0, double dx = 1.0, std::string label = {})
        : Signal(p, dtype_of<T>::value, false, count, x0, dx, std::move(label)) {}

    template <class T>
    Signal(const std::complex<T>* p, std::size_t count, double x0 = 0.0, double dx = 1.0,
           std::string label = {})
        : Signal(p, dtype_of<T>::value, true, count, x0, dx, std::move(label)) {}

    template <class T>
    Signal(const std::vector<T>& v, double x0 = 0.0, double dx = 1.0, std::string label = {})
        : Signal(v.data(), v.size(), x0, dx, std::move(label)) {}

    double x(std::size_t i) const { return xstart + static_cast<double>(i) * xdelta; }
};

}  // namespace ssp
