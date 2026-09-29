#include "ssp/render/axis.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace ssp {
namespace {

/// Fewest decimals that print `step` exactly (steps are 1, 2, 2.5 or 5 × 10^k).
int decimals_for(double step) {
    int d = std::max(0, static_cast<int>(-std::floor(std::log10(step) + 1e-9)));
    const double scaled = step * std::pow(10.0, d);
    if (std::abs(scaled - std::round(scaled)) > 1e-6 * scaled) ++d;
    return std::min(d, 15);
}

std::string fixed(double v, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", decimals, v);
    std::string s = buf;
    // "-0", "-0.00" → "0", "0.00"
    if (s[0] == '-' && s.find_first_not_of("-0.") == std::string::npos) s.erase(0, 1);
    return s;
}

std::string mult_note(double mult) {
    if (mult == 1.0) return {};
    return "x1e" + std::to_string(static_cast<int>(std::lround(std::log10(mult))));
}

}  // namespace

Tics nice_tics(double dmin, double dmax, int ndiv) {
    if (dmax == dmin) return {1.0, dmin};
    const double dran = std::abs(dmax - dmin);
    const double df = dran / std::max(1, ndiv);
    double sig = std::log10(std::max(df, 1.0e-36));
    const double nsig = sig < 0.0 ? std::ceil(sig) - 1.0 : std::floor(sig);
    const double ddf = df * std::pow(10.0, -nsig);
    sig = std::pow(10.0, nsig);

    double dtic;
    if (ddf < 1.75) dtic = sig;
    else if (ddf < 2.25) dtic = 2.0 * sig;
    else if (ddf < 3.5) dtic = 2.5 * sig;
    else if (ddf < 7.0) dtic = 5.0 * sig;
    else dtic = 10.0 * sig;
    if (dtic == 0.0) dtic = 1.0;

    Tics t;
    if (dmax >= dmin) {
        const double nseg = std::floor(dmin >= 0.0 ? dmin / dtic + 0.995 : dmin / dtic - 0.005);
        t = {dtic, nseg * dtic};
    } else {
        const double nseg = std::floor(dmin >= 0.0 ? dmin / dtic + 0.005 : dmin / dtic - 0.995);
        t = {-dtic, nseg * dtic};
    }
    if (t.dtic1 + t.dtic == t.dtic1) t.dtic = dmax - dmin;
    return t;
}

double eng_mult(double a, double b) {
    const double absmax = std::max(std::abs(a), std::abs(b));
    if (absmax == 0.0 || (absmax >= 1e-3 && absmax < 1e4)) return 1.0;
    int k = static_cast<int>(std::floor(std::log10(absmax) / 3.0));
    return std::pow(10.0, 3 * k);
}

std::size_t AxisTicks::max_label_chars() const {
    std::size_t n = 0;
    for (const auto& l : labels) n = std::max(n, l.size());
    return n;
}

AxisTicks make_ticks(double lo, double hi, int ndiv) {
    AxisTicks out;
    if (!(hi > lo) || !std::isfinite(lo) || !std::isfinite(hi)) return out;
    const Tics t = nice_tics(lo, hi, ndiv);
    const double eps = 1e-9 * (hi - lo);
    for (int k = 0; k < 1000; ++k) {
        double v = t.dtic1 + k * t.dtic;
        if (v > hi + eps) break;
        if (v < lo - eps) continue;
        if (std::abs(v) < 1e-9 * t.dtic) v = 0.0;
        out.values.push_back(v);
    }
    if (out.values.empty()) return out;

    auto label_all = [&](double offset, double mult) {
        const int d = decimals_for(t.dtic / mult);
        out.labels.clear();
        for (double v : out.values) out.labels.push_back(fixed((v - offset) / mult, d));
    };

    const double mult = eng_mult(lo, hi);
    label_all(0.0, mult);
    out.note = mult_note(mult);
    if (out.max_label_chars() <= static_cast<std::size_t>(kMaxTickLabelChars)) return out;

    // Too many digits to tell ticks apart (e.g. 1e9 + small span): label relative to the
    // first tick and show that offset once.
    const double offset = out.values.front();
    const double rmult = eng_mult(0.0, hi - lo);
    label_all(offset, rmult);
    out.note = "+" + format_g(offset, 12);
    if (rmult != 1.0) out.note += " " + mult_note(rmult);
    return out;
}

std::string format_g(double v, int sig) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", sig, v);
    return buf;
}

}  // namespace ssp
