#include "ssp/plot.h"

#include "render/png.h"
#include "render/trace_plot.h"

namespace ssp {

void save_png(const std::vector<Signal>& signals, const PlotOptions& options,
              const std::string& path, int width, int height) {
    if (width < 64 || height < 64) throw std::invalid_argument("save_png: image too small");
    TracePlot plot(signals, options);
    Framebuffer fb(width, height);
    plot.render(fb);
    write_png(fb, path);
}

}  // namespace ssp
