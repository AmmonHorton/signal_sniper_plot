/// @file scan.h
/// @brief The one tight loop over raw samples, specialised per element type and component.
#pragma once

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include "core/cancel.h"
#include "core/component.h"
#include "core/dispatch.h"
#include "ssp/types.h"

namespace ssp {

constexpr std::size_t kCancelCheckEvery = 1 << 16;  // samples between cancellation checks

namespace detail {

template <class T, bool Cplx, Comp C, class F>
void for_each_t(const T* p, std::ptrdiff_t stride, std::size_t i0, std::size_t i1,
                const CancelToken& ct, F& f) {
    const std::ptrdiff_t step = Cplx ? 2 * stride : stride;
    for (std::size_t chunk = i0; chunk < i1; chunk += kCancelCheckEvery) {
        ct.check();
        const std::size_t end = std::min(i1, chunk + kCancelCheckEvery);
        const T* q = p + static_cast<std::ptrdiff_t>(chunk) * step;
        for (std::size_t i = chunk; i < end; ++i, q += step) {
            const double re = static_cast<double>(q[0]);
            const double im = Cplx ? static_cast<double>(q[1]) : 0.0;
            f(comp_value<C>(re, im));
        }
    }
}

template <class T, bool Cplx, class F>
void for_each_pair_t(const T* p, std::ptrdiff_t stride, std::size_t i0, std::size_t i1,
                     const CancelToken& ct, F& f) {
    const std::ptrdiff_t step = Cplx ? 2 * stride : stride;
    for (std::size_t chunk = i0; chunk < i1; chunk += kCancelCheckEvery) {
        ct.check();
        const std::size_t end = std::min(i1, chunk + kCancelCheckEvery);
        const T* q = p + static_cast<std::ptrdiff_t>(chunk) * step;
        for (std::size_t i = chunk; i < end; ++i, q += step) {
            f(static_cast<double>(q[0]), Cplx ? static_cast<double>(q[1]) : 0.0);
        }
    }
}

}  // namespace detail

/// @brief Call f(value) with component c of every sample in [i0, i1), in order.
/// @throws Cancelled when `ct` is cancelled.
template <class F>
void for_each_value(const Signal& s, Comp c, std::size_t i0, std::size_t i1,
                    const CancelToken& ct, F&& f) {
    i1 = std::min(i1, s.n);
    if (i0 >= i1) return;
    dispatch(s.dtype, [&](auto tag) {
        using T = typename decltype(tag)::type;
        const T* p = static_cast<const T*>(s.data);
        auto run = [&](auto cplx) {
            constexpr bool X = decltype(cplx)::value;
            switch (c) {
                case Comp::Re:    return detail::for_each_t<T, X, Comp::Re>(p, s.stride, i0, i1, ct, f);
                case Comp::Im:    return detail::for_each_t<T, X, Comp::Im>(p, s.stride, i0, i1, ct, f);
                case Comp::Mag:   return detail::for_each_t<T, X, Comp::Mag>(p, s.stride, i0, i1, ct, f);
                case Comp::Phase: return detail::for_each_t<T, X, Comp::Phase>(p, s.stride, i0, i1, ct, f);
            }
        };
        if (s.complex) {
            run(std::true_type{});
        } else {
            run(std::false_type{});
        }
    });
}

/// @brief Call f(re, im) for every sample in [i0, i1), in order (im = 0 for real data).
/// @throws Cancelled when `ct` is cancelled.
template <class F>
void for_each_sample(const Signal& s, std::size_t i0, std::size_t i1, const CancelToken& ct,
                     F&& f) {
    i1 = std::min(i1, s.n);
    if (i0 >= i1) return;
    dispatch(s.dtype, [&](auto tag) {
        using T = typename decltype(tag)::type;
        const T* p = static_cast<const T*>(s.data);
        if (s.complex) {
            detail::for_each_pair_t<T, true>(p, s.stride, i0, i1, ct, f);
        } else {
            detail::for_each_pair_t<T, false>(p, s.stride, i0, i1, ct, f);
        }
    });
}

}  // namespace ssp
