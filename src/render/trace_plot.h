/// @file trace_plot.h
/// @brief xplot scene: signals + options → autoscale → reduce → paint.
#pragma once

#include <memory>
#include <vector>

#include "core/lod.h"
#include "render/frame.h"
#include "render/trace.h"
#include "ssp/plot.h"

namespace ssp {

class TracePlot {
public:
    /// @throws std::invalid_argument for unusable signals or options.
    TracePlot(std::vector<Signal> signals, PlotOptions options);

    /// @brief Full synchronous render: autoscale, reduce every visible trace, paint.
    void render(Framebuffer& fb, const CancelToken& ct = {});

    CMode cmode() const { return cmode_; }
    Range xview() const { return xview_; }
    Range yview() const { return yview_; }

    // The steps render() is made of, exposed for the interactive (phase 2) pipeline.
    Range autoscale_x() const;
    void reduce(const XView& view, const CancelToken& ct = {});
    Range autoscale_y() const;
    void paint(Framebuffer& fb, const XView& view, Range y);

private:
    double xstart(std::size_t i) const { return opts_.index ? 0.0 : sigs_[i].xstart; }
    double xdelta(std::size_t i) const { return opts_.index ? 1.0 : sigs_[i].xdelta; }
    XView make_view(Range x, int width) const;

    std::vector<Signal> sigs_;
    PlotOptions opts_;
    CMode cmode_;
    bool show_legend_ = false;  // only when there is more than one trace or a trace has a name
    std::vector<std::unique_ptr<Lod>> lods_;
    std::vector<TraceBins> bins_;
    Theme theme_;
    Range xview_, yview_;
};

}  // namespace ssp
