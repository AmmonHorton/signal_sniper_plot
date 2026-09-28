/// @file component.h
/// @brief Signal components (re, im, |x|, arg x) and the display transform for each CMode.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "ssp/types.h"

namespace ssp {

/// @brief The quantity a min/max pyramid is built over. Log modes reuse Mag because
/// log10 is monotonic, so min/max commute with it.
enum class Comp : uint8_t { Re, Im, Mag, Phase };
constexpr int kNumComps = 4;

inline Comp comp_of(CMode m) {
    switch (m) {
        case CMode::Real:  return Comp::Re;
        case CMode::Imag:  return Comp::Im;
        case CMode::Phase: return Comp::Phase;
        case CMode::Mag:
        case CMode::Log10:
        case CMode::Log20: return Comp::Mag;
        case CMode::Auto:
        case CMode::IR:    break;
    }
    throw std::invalid_argument("cmode has no single y component");
}

template <Comp C>
inline double comp_value(double re, double im) {
    if constexpr (C == Comp::Re) return re;
    if constexpr (C == Comp::Im) return im;
    if constexpr (C == Comp::Mag) return std::sqrt(re * re + im * im);
    if constexpr (C == Comp::Phase) return std::atan2(im, re);
}

/// @brief Magnitudes below this plot as the floor instead of -inf in log modes.
constexpr double kLogFloor = 1e-20;

/// @brief Map a component value to plotted y. Monotonic non-decreasing for every mode, so
/// it can be applied to the ends of a min/max span.
inline double display_value(CMode m, PhaseUnits units, double v) {
    switch (m) {
        case CMode::Log10: return 10.0 * std::log10(std::max(v, kLogFloor));
        case CMode::Log20: return 20.0 * std::log10(std::max(v, kLogFloor));
        case CMode::Phase:
            if (units == PhaseUnits::Degrees) return v * (180.0 / M_PI);
            if (units == PhaseUnits::Cycles) return v / (2.0 * M_PI);
            return v;
        default: return v;
    }
}

/// @brief Running min/max of finite values. Empty until the first finite value.
struct Span {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();

    bool empty() const { return !(lo <= hi); }
    void add(double v) {
        if (!std::isfinite(v)) return;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
    }
    void merge(const Span& o) {
        if (o.lo < lo) lo = o.lo;
        if (o.hi > hi) hi = o.hi;
    }
};

}  // namespace ssp
