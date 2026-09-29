#include <memory>
#include <stdexcept>

#include "ssp/app/compose.h"
#include "ssp/app/raster_content.h"
#include "ssp/render/png.h"
#include "ssp/plot.h"

namespace ssp {

void save_png(const std::vector<Signal>& signals, const PlotOptions& options,
              const std::string& path, int width, int height) {
    if (width < 64 || height < 64) throw std::invalid_argument("save_png: image too small");
    Screen s(std::make_unique<TraceContent>(signals, options));
    Framebuffer fb(width, height);
    render_headless(fb, s);
    write_png(fb, path);
}

void save_raster_png(const Signal& data, const RasterOptions& options, const std::string& path,
                     int width, int height) {
    if (width < 64 || height < 64) throw std::invalid_argument("save_raster_png: image too small");
    Screen s(std::make_unique<RasterContent>(data, options));
    Framebuffer fb(width, height);
    render_headless(fb, s);
    write_png(fb, path);
}

}  // namespace ssp
