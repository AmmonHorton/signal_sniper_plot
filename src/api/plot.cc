#include "ssp/plot.h"

#include "app/session.h"
#include "app/trace_content.h"
#include "platform/backend.h"

namespace ssp {

Session::Session(std::vector<Signal> signals, PlotOptions options)
    : signals_(std::move(signals)), options_(std::move(options)) {
    TraceContent check(signals_, options_);  // validate now, not when the window opens
}

void Session::run() {
    auto backend = open_x11(options_.width, options_.height,
                            options_.title.empty() ? "signal_sniper_plot" : options_.title);
    run_plot(signals_, options_, *backend, interrupt_);
}

void plot(const std::vector<Signal>& signals, const PlotOptions& options) {
    Session(signals, options).run();
}

}  // namespace ssp
