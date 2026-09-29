#include "ssp/plot.h"

#include "ssp/app/raster_content.h"
#include "ssp/app/session.h"
#include "ssp/app/trace_content.h"
#include "ssp/platform/backend.h"

namespace ssp {

Session::Session(std::vector<Signal> signals, PlotOptions options)
    : signals_(std::move(signals)), options_(std::move(options)) {
    TraceContent check(signals_, options_);  // validate now, not when the window opens
}

Session::Session(Signal data, RasterOptions options) : signals_{std::move(data)} {
    RasterContent check(signals_[0], options);
    options_.title = options.title;
    options_.width = options.width;
    options_.height = options.height;
    raster_ = std::move(options);
}

void Session::run() {
    auto backend = open_x11(options_.width, options_.height,
                            options_.title.empty() ? "signal_sniper_plot" : options_.title);
    std::unique_ptr<Content> content;
    if (raster_) {
        content = std::make_unique<RasterContent>(signals_[0], *raster_);
    } else {
        content = std::make_unique<TraceContent>(signals_, options_);
    }
    run_session(std::move(content), *backend, interrupt_);
}

void plot(const std::vector<Signal>& signals, const PlotOptions& options) {
    Session(signals, options).run();
}

void raster(const Signal& data, const RasterOptions& options) { Session(data, options).run(); }

}  // namespace ssp
