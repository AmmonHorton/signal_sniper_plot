/// @file cancel.h
/// @brief Cooperative cancellation for long-running reduce work.
#pragma once

#include <atomic>
#include <cstdint>

namespace ssp {

/// @brief Thrown by CancelToken::check() once the job it belongs to is superseded.
struct Cancelled {};

/// @brief A job is cancelled when the shared generation counter moves past the job's own
/// generation. A default-constructed token never cancels.
class CancelToken {
public:
    CancelToken() = default;
    CancelToken(const std::atomic<uint64_t>& generation, uint64_t mine)
        : gen_(&generation), mine_(mine) {}

    bool cancelled() const {
        return gen_ && gen_->load(std::memory_order_relaxed) != mine_;
    }
    void check() const {
        if (cancelled()) throw Cancelled{};
    }

private:
    const std::atomic<uint64_t>* gen_ = nullptr;
    uint64_t mine_ = 0;
};

}  // namespace ssp
