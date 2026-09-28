/// @file lod.h
/// @brief Level-of-detail min/max pyramid so any 1-D view costs O(pixels), not O(samples).
#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "core/cancel.h"
#include "core/component.h"
#include "ssp/types.h"

namespace ssp {

/// @brief Component value of sample i.
double sample_value(const Signal& s, Comp c, std::size_t i);

/// @brief Min/max of component c over samples [i0, i1) by reading every sample.
Span scan(const Signal& s, Comp c, std::size_t i0, std::size_t i1,
          const CancelToken& ct = {});

/// @brief Answers "min/max of component c over samples [i0, i1)" for one Signal.
///
/// Level 0 holds one span per block of kBlock samples; each level above merges kFan
/// children. Level-0 blocks are filled as a side effect of whatever queries scan them
/// (tracked per block), so the first full-extent view builds the pyramid and every later
/// view reuses it. Upper levels are built once level 0 is complete.
///
/// Not thread-safe: use from one thread (the render worker).
class Lod {
public:
    static constexpr std::size_t kBlock = 1024;
    static constexpr std::size_t kFan = 32;

    explicit Lod(Signal s) : sig_(std::move(s)) {}

    Span query(Comp c, std::size_t i0, std::size_t i1, const CancelToken& ct = {});

    /// @brief Fraction of level-0 blocks built for component c, in [0, 1].
    double coverage(Comp c) const;

    const Signal& signal() const { return sig_; }

private:
    struct Pyramid {
        std::vector<std::vector<Span>> levels;  // levels[0] = per-block spans
        std::vector<uint8_t> valid;             // per level-0 block
        std::size_t nvalid = 0;
    };

    Pyramid& pyramid(Comp c);
    Span query_blocks(Pyramid& p, Comp c, std::size_t b0, std::size_t b1, const CancelToken& ct);
    void build_upper(Pyramid& p);

    Signal sig_;
    std::array<std::unique_ptr<Pyramid>, kNumComps> pyr_;
};

}  // namespace ssp
