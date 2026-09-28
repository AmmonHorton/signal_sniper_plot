#include <memory>
#include <stdexcept>

#include "app/compose.h"
#include "render/png.h"
#include "ssp/plot.h"

namespace ssp {

void save_png(const std::vector<Signal>& signals, const PlotOptions& options,
              const std::string& path, int width, int height) {
    if (width < 64 || height < 64) throw std::invalid_argument("save_png: image too small");
    Screen s(std::make_unique<TraceContent>(signals, options), options);
    Framebuffer fb(width, height);
    render_headless(fb, s);
    write_png(fb, path);
}

}  // namespace ssp
