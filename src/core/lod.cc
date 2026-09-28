#include "core/lod.h"

#include <algorithm>
#include <cmath>

#include "core/dispatch.h"
#include "core/scan.h"

namespace ssp {
namespace {

/// For real data, Im and Phase follow exactly from the real span (Im is 0; the phase is 0
/// or pi depending on sign), so they need no pass over the samples. |x| does not: a span
/// of [-5, 3] says nothing about how close to 0 the samples get, so Mag has its own pyramid.
bool derivable_from_real(Comp c) { return c == Comp::Im || c == Comp::Phase; }

Span derive_from_real(Comp c, const Span& re) {
    Span out;
    if (re.empty()) return out;
    if (c == Comp::Im) {
        out.add(0.0);
    } else if (c == Comp::Phase) {
        if (re.hi >= 0.0) out.add(0.0);
        if (re.lo < 0.0) out.add(M_PI);
    }
    return out;
}

}  // namespace

double sample_value(const Signal& s, Comp c, std::size_t i) {
    return dispatch(s.dtype, [&](auto tag) {
        using T = typename decltype(tag)::type;
        const T* p = static_cast<const T*>(s.data);
        const std::ptrdiff_t off = static_cast<std::ptrdiff_t>(i) * s.stride * (s.complex ? 2 : 1);
        const double re = static_cast<double>(p[off]);
        const double im = s.complex ? static_cast<double>(p[off + 1]) : 0.0;
        switch (c) {
            case Comp::Re:    return comp_value<Comp::Re>(re, im);
            case Comp::Im:    return comp_value<Comp::Im>(re, im);
            case Comp::Mag:   return comp_value<Comp::Mag>(re, im);
            case Comp::Phase: return comp_value<Comp::Phase>(re, im);
        }
        return re;
    });
}

Span scan(const Signal& s, Comp c, std::size_t i0, std::size_t i1, const CancelToken& ct) {
    Span out;
    for_each_value(s, c, i0, i1, ct, [&out](double v) { out.add(v); });
    return out;
}

Lod::Pyramid& Lod::pyramid(Comp c) {
    auto& slot = pyr_[static_cast<int>(c)];
    if (!slot) {
        // Sized for the current length. Live data will need to grow the last level-0
        // block and the tail of each level; nothing here assumes n is final beyond this.
        const std::size_t nblocks = (sig_.n + kBlock - 1) / kBlock;
        slot = std::make_unique<Pyramid>();
        slot->levels.emplace_back(nblocks);
        slot->valid.assign(nblocks, 0);
    }
    return *slot;
}

double Lod::coverage(Comp c) const {
    const auto& slot = pyr_[static_cast<int>(c)];
    if (!slot || slot->valid.empty()) return 0.0;
    return static_cast<double>(slot->nvalid) / static_cast<double>(slot->valid.size());
}

void Lod::build_upper(Pyramid& p) {
    while (p.levels.back().size() > 1) {
        const auto& below = p.levels.back();
        std::vector<Span> up((below.size() + kFan - 1) / kFan);
        for (std::size_t i = 0; i < below.size(); ++i) up[i / kFan].merge(below[i]);
        p.levels.push_back(std::move(up));
    }
}

Span Lod::query_blocks(Pyramid& p, Comp c, std::size_t b0, std::size_t b1,
                       const CancelToken& ct) {
    Span out;
    auto& level0 = p.levels[0];

    if (p.nvalid < level0.size()) {
        // Pyramid incomplete: use built blocks, scan (and store) the rest.
        for (std::size_t b = b0; b < b1; ++b) {
            if (!p.valid[b]) {
                level0[b] = scan(sig_, c, b * kBlock, (b + 1) * kBlock, ct);
                p.valid[b] = 1;
                ++p.nvalid;
            }
            out.merge(level0[b]);
        }
        if (p.nvalid == level0.size()) build_upper(p);
        return out;
    }

    // Complete: merge the ragged ends at each level and climb for the aligned middle.
    std::size_t a = b0, b = b1;
    for (std::size_t k = 0; a < b; ++k) {
        const auto& lv = p.levels[k];
        if (k + 1 == p.levels.size()) {
            for (std::size_t i = a; i < b; ++i) out.merge(lv[i]);
            break;
        }
        const std::size_t a2 = (a + kFan - 1) / kFan;
        const std::size_t b2 = (b == lv.size()) ? p.levels[k + 1].size() : b / kFan;
        if (a2 >= b2) {
            for (std::size_t i = a; i < b; ++i) out.merge(lv[i]);
            break;
        }
        for (std::size_t i = a; i < a2 * kFan; ++i) out.merge(lv[i]);
        for (std::size_t i = b2 * kFan; i < b; ++i) out.merge(lv[i]);
        a = a2;
        b = b2;
    }
    return out;
}

Span Lod::query(Comp c, std::size_t i0, std::size_t i1, const CancelToken& ct) {
    i1 = std::min(i1, sig_.n);
    if (i0 >= i1) return {};
    if (!sig_.complex && derivable_from_real(c)) return derive_from_real(c, query(Comp::Re, i0, i1, ct));

    // Blocks fully inside [i0, i1); the last block may be short.
    const std::size_t nblocks = (sig_.n + kBlock - 1) / kBlock;
    const std::size_t b0 = (i0 + kBlock - 1) / kBlock;
    const std::size_t b1 = (i1 == sig_.n) ? nblocks : i1 / kBlock;
    if (b0 >= b1) return scan(sig_, c, i0, i1, ct);

    Span out = scan(sig_, c, i0, b0 * kBlock, ct);
    out.merge(query_blocks(pyramid(c), c, b0, b1, ct));
    out.merge(scan(sig_, c, std::min(b1 * kBlock, i1), i1, ct));
    return out;
}

}  // namespace ssp
