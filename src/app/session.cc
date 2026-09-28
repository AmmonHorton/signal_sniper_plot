#include "app/session.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <ctime>
#include <memory>
#include <string>

#include "app/compose.h"
#include "app/controller.h"
#include "app/job_runner.h"
#include "app/overlay.h"
#include "app/screen.h"
#include "render/png.h"

namespace ssp {
namespace {

using Clock = std::chrono::steady_clock;

constexpr int kFrameMs = 16;       // redraw rate while a reduce is streaming in
constexpr int kInterruptMs = 100;  // how often `interrupt` is polled when idle
/// After a view change, keep showing the old plot this long before showing an empty one;
/// most reduces finish well within it once the pyramid exists.
constexpr auto kBlankDelay = std::chrono::milliseconds(60);

std::string png_name(const std::string& title) {
    std::string stem;
    for (char c : title) stem += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
    char ts[32];
    const std::time_t now = std::time(nullptr);
    std::strftime(ts, sizeof ts, "%Y%m%d-%H%M%S", std::localtime(&now));
    return "ssp_" + (stem.empty() ? std::string("plot") : stem) + "_" + ts + ".png";
}

void copy_rect(const Framebuffer& from, Framebuffer& to, Rect r) {
    r = r.intersect(from.bounds());
    for (int y = r.y; y < r.bottom(); ++y) {
        std::copy(from.row(y) + r.x, from.row(y) + r.right(), to.row(y) + r.x);
    }
}

class Loop {
public:
    Loop(std::vector<Signal> signals, const PlotOptions& o, Backend& be)
        : be_(be),
          screen_(std::make_unique<TraceContent>(std::move(signals), o), o),
          plot_(be.width(), be.height()),
          frame_(plot_) {}

    void run(const std::function<bool()>& interrupt) {
        submit();
        std::vector<InputEvent> events;
        for (;;) {
            // Never block while finished columns are waiting to be shown.
            const int timeout = done() != shown_done_ && done() > 0 ? 0
                                : rendering()                     ? kFrameMs
                                : interrupt                       ? kInterruptMs
                                                                  : -1;
            events.clear();
            be_.wait_events(events, timeout, jobs_.wake_fd());
            jobs_.drain();

            unsigned inv = kNone;
            std::vector<Rect> damage;
            for (const InputEvent& e : events) {
                if (e.type == InputEvent::Type::Resize) {
                    plot_ = Framebuffer(e.w, e.h);
                    frame_ = plot_;
                    first_ = true;  // show the new layout right away
                    inv |= kReduce;
                } else if (e.type == InputEvent::Type::Exposed) {
                    damage.push_back(e.rect);
                } else {
                    inv |= handle_event(screen_, e, rendering());
                    if (inv & kCancel) {
                        stop();
                        inv &= ~kCancel;
                    }
                }
            }
            if ((inv & kQuit) || (interrupt && interrupt())) break;
            if (inv & (kReduce | kResume)) submit();
            if (inv & kSave) save();
            frame(inv, damage);
        }
        jobs_.cancel();
        jobs_.wait_idle();
    }

private:
    bool rendering() const {
        return screen_.result && !screen_.stopped &&
               !screen_.result->ended.load(std::memory_order_acquire);
    }
    int done() const {
        return screen_.result ? screen_.result->done.load(std::memory_order_acquire) : 0;
    }
    void submit() {
        Screen& s = screen_;
        s.layout = compute_layout(plot_.width(), plot_.height());
        const XView view = s.xview();
        View& level = s.views.top_mut();
        if (level.cache && s.content->reusable(*level.cache, view, s.set.index, level.y)) {
            jobs_.cancel();  // e.g. unzoom while the deeper level was still rendering
            s.result = level.cache;
            s.stopped = false;
            shown_done_ = -1;  // repaint from it now
            return;
        }
        auto r = s.content->new_result(view, s.set.index, level.y);
        level.cache = r;
        s.result = r;
        s.stopped = false;
        shown_done_ = 0;  // keep the old picture until columns arrive (or kBlankDelay passes)
        submitted_ = Clock::now();
        const TraceContent* content = s.content.get();
        JobRunner* jobs = &jobs_;
        jobs_.submit([content, r, jobs](const CancelToken& ct) {
            struct Ended {
                ReduceResult& r;
                ~Ended() { r.ended.store(true, std::memory_order_release); }
            } ended{*r};
            content->reduce(*r, ct, [jobs] { jobs->notify(); });
        });
    }

    void stop() {
        jobs_.cancel();
        screen_.stopped = true;
    }

    void save() {
        const std::string path = png_name(screen_.set.title);
        try {
            write_png(plot_, path);
            screen_.message = "saved " + path;
        } catch (const std::exception& e) {
            screen_.message = e.what();
        }
    }

    /// Recompose if the data changed, redraw overlays, and upload what changed.
    void frame(unsigned inv, std::vector<Rect>& damage) {
        const int d = done();
        const bool new_columns = d != shown_done_;
        const bool blank_due = shown_done_ == 0 && d == 0 && !shown_blank_ &&
                               Clock::now() - submitted_ >= kBlankDelay;
        if ((inv & kRepaint) || new_columns || blank_due || first_) {
            compose(plot_, screen_, d, theme_);
            shown_done_ = d;
            shown_blank_ = d == 0;
            first_ = false;
            frame_ = plot_;
            overlay_.clear();
            damage.assign(1, plot_.bounds());
        } else {
            for (const Rect& r : overlay_) copy_rect(plot_, frame_, r);
            damage.insert(damage.end(), overlay_.begin(), overlay_.end());
        }
        if (d > 0) shown_blank_ = false;

        const double progress =
            screen_.result ? static_cast<double>(d) / screen_.result->width() : 1.0;
        overlay_ = draw_overlay(frame_, screen_, theme_, rendering(), progress);
        damage.insert(damage.end(), overlay_.begin(), overlay_.end());
        be_.present(frame_, damage);
    }

    Backend& be_;
    Theme theme_;
    Screen screen_;
    Framebuffer plot_;   // composed plot, no overlays (also what Ctrl-S saves)
    Framebuffer frame_;  // plot_ + overlays: what is on screen
    std::vector<Rect> overlay_;
    JobRunner jobs_;     // declared last: destroyed (joined) first
    int shown_done_ = -1;
    bool shown_blank_ = false;
    bool first_ = true;
    Clock::time_point submitted_;
};

}  // namespace

void run_plot(std::vector<Signal> signals, const PlotOptions& options, Backend& backend,
              const std::function<bool()>& interrupt) {
    Loop(std::move(signals), options, backend).run(interrupt);
}

}  // namespace ssp
