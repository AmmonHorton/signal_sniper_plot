/// @file axis.h
/// @brief "Nice" tick placement (port of SigPlot mx.tics) and tick/readout number formatting.
#pragma once

#include <string>
#include <vector>

namespace ssp {

/// @brief First tick at or below dmin and the tick spacing, splitting [dmin, dmax] into
/// about `ndiv` steps of 1, 2, 2.5 or 5 × 10^k. Port of SigPlot `mx.tics` (non-timecode).
struct Tics {
    double dtic = 1.0;   ///< Spacing.
    double dtic1 = 0.0;  ///< First tick.
};
Tics nice_tics(double dmin, double dmax, int ndiv);

/// @brief Engineering multiplier (10^(3k)) for labelling [a, b]; 1 for ordinary magnitudes.
double eng_mult(double a, double b);

/// @brief Tick labels are kept to this many characters (longer spans switch to offset labels).
constexpr int kMaxTickLabelChars = 10;

/// @brief Ticks inside [lo, hi] with labels short enough to stay unique.
struct AxisTicks {
    std::vector<double> values;
    std::vector<std::string> labels;
    /// Shown once near the axis when labels are scaled or offset, e.g. "x1e-6" or "+1.23e9 x1e3".
    std::string note;
    std::size_t max_label_chars() const;
};
AxisTicks make_ticks(double lo, double hi, int ndiv);

/// @brief `%.{sig}g`, for readouts.
std::string format_g(double v, int sig = 9);

}  // namespace ssp
