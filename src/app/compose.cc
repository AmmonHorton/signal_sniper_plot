#include "app/compose.h"

#include "render/colormap.h"

namespace ssp {

void compose(Framebuffer& fb, Screen& s, int done, const Theme& th, bool mode_label) {
    s.layout = s.layout_for(fb.width(), fb.height());
    s.xshown = s.x_for(done);
    s.yshown = s.y_for(done);
    const Range x = s.xshown;
    const Range y = s.yshown;
    const Rect& p = s.layout.plot;

    fb.fill(th.bg);
    const AxisTicks xt = make_ticks(x.lo, x.hi, xdivisions(p.w));
    const AxisTicks yt = make_ticks(y.lo, y.hi, ydivisions(p.h));
    const bool raster = s.content->is_raster();
    // Rasters cover the grid, so theirs is drawn again on top of the image below.
    draw_axes(fb, s.layout, x.lo, x.hi, xt, y.lo, y.hi, yt, s.set.grid && !raster, th, s.ydown());
    s.zshown = s.z_for(done);
    if (s.result) s.content->paint(fb, *s.result, PaintArgs{p, done, x, y, s.zshown, s.set, th});
    if (raster) {
        if (s.set.grid) draw_axes(fb, s.layout, x.lo, x.hi, xt, y.lo, y.hi, yt, true, th, true);
        draw_colorbar(fb, s.layout, s.zshown.lo, s.zshown.hi, colormap_lut(s.set.cmap), th);
    }
    fb.rect_outline({p.x - 1, p.y - 1, p.w + 2, p.h + 2}, th.fg);

    s.legend_hits.clear();
    if (s.set.legend) s.legend_hits = draw_legend(fb, p, s.content->legend(th), th);
    draw_title(fb, s.layout, s.set.title, th);
    if (mode_label) fb.text(Framebuffer::kCharW, s.layout.readout.y + 4, cmode_name(s.set.cmode), th.fg);
}

void render_headless(Framebuffer& fb, Screen& s, const Theme& th) {
    s.layout = s.layout_for(fb.width(), fb.height());
    s.result = s.content->new_result(s.request());
    s.content->reduce(*s.result, {});
    s.result->ended = true;
    compose(fb, s, s.result->total, th, /*mode_label=*/true);
}

}  // namespace ssp
