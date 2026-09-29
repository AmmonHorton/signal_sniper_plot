#include "ssp/app/job_runner.h"

#include <fcntl.h>
#include <unistd.h>

#include <cstdio>
#include <stdexcept>

namespace ssp {

JobRunner::JobRunner() {
    if (::pipe(pipe_) != 0) throw std::runtime_error("JobRunner: pipe() failed");
    for (int fd : pipe_) ::fcntl(fd, F_SETFL, ::fcntl(fd, F_GETFL) | O_NONBLOCK);
    thread_ = std::thread([this] { loop(); });
}

JobRunner::~JobRunner() {
    {
        std::lock_guard<std::mutex> lk(mu_);
        stop_ = true;
        ++gen_;
    }
    cv_.notify_all();
    thread_.join();
    ::close(pipe_[0]);
    ::close(pipe_[1]);
}

void JobRunner::submit(Job job) {
    {
        std::lock_guard<std::mutex> lk(mu_);
        pending_ = std::move(job);
        pending_gen_ = ++gen_;  // also cancels the running job
    }
    cv_.notify_all();
}

void JobRunner::cancel() {
    std::lock_guard<std::mutex> lk(mu_);
    pending_ = nullptr;
    ++gen_;
}

void JobRunner::wait_idle() {
    std::unique_lock<std::mutex> lk(mu_);
    cv_.wait(lk, [&] { return !running_ && !pending_; });
}

void JobRunner::notify() {
    if (wake_pending_.exchange(true)) return;  // one byte in flight is enough
    const char b = 1;
    (void)!::write(pipe_[1], &b, 1);
}

void JobRunner::drain() {
    char buf[64];
    wake_pending_ = false;
    while (::read(pipe_[0], buf, sizeof buf) > 0) {
    }
}

void JobRunner::loop() {
    for (;;) {
        Job job;
        uint64_t mine;
        {
            std::unique_lock<std::mutex> lk(mu_);
            cv_.wait(lk, [&] { return stop_ || pending_; });
            if (stop_) return;
            job = std::move(pending_);
            pending_ = nullptr;
            mine = pending_gen_;
            running_ = true;
        }
        try {
            job(CancelToken(gen_, mine));
        } catch (const Cancelled&) {
        } catch (const std::exception& e) {
            std::fprintf(stderr, "ssp: render job failed: %s\n", e.what());
        }
        job = nullptr;  // release captures before reporting idle
        {
            std::lock_guard<std::mutex> lk(mu_);
            running_ = false;
        }
        cv_.notify_all();
        notify();
    }
}

}  // namespace ssp
