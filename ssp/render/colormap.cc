#include "ssp/render/colormap.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ssp {
namespace {

struct Stop {
    float pos, r, g, b;
};
struct Table {
    const char* name;
    std::vector<Stop> stops;
};

const Table kTables[kNumColormaps] = {
#include "ssp/render/colormaps.inc"
};

ColorLut build(const Table& t) {
    ColorLut lut{};
    for (int i = 0; i < 256; ++i) {
        const float p = i / 255.0f;
        std::size_t k = 1;
        while (k + 1 < t.stops.size() && t.stops[k].pos < p) ++k;
        const Stop& a = t.stops[k - 1];
        const Stop& b = t.stops[k];
        const float f = b.pos > a.pos ? std::clamp((p - a.pos) / (b.pos - a.pos), 0.0f, 1.0f) : 1.0f;
        auto ch = [&](float x, float y) { return static_cast<uint32_t>(std::lround(255.0f * (x + (y - x) * f))); };
        lut[i] = ch(a.r, b.r) << 16 | ch(a.g, b.g) << 8 | ch(a.b, b.b);
    }
    return lut;
}

}  // namespace

const ColorLut& colormap_lut(Colormap c) {
    static const std::array<ColorLut, kNumColormaps> luts = [] {
        std::array<ColorLut, kNumColormaps> out{};
        for (int i = 0; i < kNumColormaps; ++i) out[i] = build(kTables[i]);
        return out;
    }();
    return luts[static_cast<int>(c)];
}

const char* colormap_name(Colormap c) { return kTables[static_cast<int>(c)].name; }

}  // namespace ssp
