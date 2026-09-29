/// @file job_runner.h
/// @brief One background thread that runs the latest submitted job; newer jobs cancel older ones.
#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>

#include "ssp/core/cancel.h"

namespace ssp {

class JobRunner {
public:
    using Job = std::function<void(const CancelToken&)>;

    JobRunner();
    ~JobRunner();  ///< Cancels and joins.
    JobRunner(const JobRunner&) = delete;
    JobRunner& operator=(const JobRunner&) = delete;

    /// @brief Cancel whatever is running or queued and queue `job`. A job that throws
    /// Cancelled simply ends; `on_exit` in the job itself is the caller's business.
    void submit(Job job);
    void cancel();
    /// @brief Block until nothing is running or queued.
    void wait_idle();

    /// @brief Readable (for poll) after notify(); drain with drain().
    int wake_fd() const { return pipe_[0]; }
    /// @brief Wake the UI thread. Safe from any thread; cheap to call often.
    void notify();
    void drain();

private:
    void loop();

    std::mutex mu_;
    std::condition_variable cv_;
    Job pending_;
    uint64_t pending_gen_ = 0;
    bool running_ = false;
    bool stop_ = false;
    std::atomic<uint64_t> gen_{0};
    std::atomic<bool> wake_pending_{false};
    int pipe_[2] = {-1, -1};
    std::thread thread_;
};

}  // namespace ssp
