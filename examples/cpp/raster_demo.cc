// Interactive xraster demo: a synthetic spectrogram (noise floor, a wandering chirp and a
// steady tone), shown in 20*log10. Try x / y over the chirp, Esc to return, c for colormaps.
//   bazel run -c opt //examples/cpp:raster_demo -- [frames] [bins]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "ssp/plot.h"

int main(int argc, char** argv) {
    const std::size_t frames = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 20'000;
    const std::size_t bins = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 4096;
    std::vector<float> spec(frames * bins);
    uint32_t st = 1;
    for (std::size_t r = 0; r < frames; ++r) {
        const double t = double(r) / frames;
        const double chirp = bins * (0.1 + 0.6 * t + 0.1 * std::sin(40 * t));
        for (std::size_t b = 0; b < bins; ++b) {
            st = st * 1664525u + 1013904223u;
            const double d1 = (b - chirp) / 4.0, d2 = (b - 0.85 * bins) / 1.5;
            spec[r * bins + b] = float(1e-3 * (0.2 + (st >> 8) / double(1u << 24)) +
                                       std::exp(-d1 * d1) + 0.1 * std::exp(-d2 * d2));
        }
    }
    ssp::RasterOptions o;
    o.title = "ssp raster demo - " + std::to_string(frames) + " x " + std::to_string(bins);
    o.subsize = bins;
    o.cmode = ssp::CMode::Log20;
    o.ydelta = 1e-3;  // 1 ms per frame
    try {
        ssp::raster(ssp::Signal(spec, -0.5 * bins * 1e3, 1e3), o);  // bins 1 kHz apart, centred
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ssp_raster_demo: %s\n", e.what());
        return 1;
    }
}
