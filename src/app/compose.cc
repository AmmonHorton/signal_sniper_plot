#include "app/compose.h"

namespace ssp {

void compose(Framebuffer& fb, Screen& s, int done, const Theme& th, bool mode_label) {
    s.layout = compute_layout(fb.width(), fb.height());
    s.yshown = s.y_for(done);
    const Range x = s.views.top().x;
    const Range y = s.yshown;
    const Rect& p = s.layout.plot;

    fb.fill(th.bg);
    const AxisTicks xt = make_ticks(x.lo, x.hi, xdivisions(p.w));
    const AxisTicks yt = make_ticks(y.lo, y.hi, ydivisions(p.h));
    draw_axes(fb, s.layout, x.lo, x.hi, xt, y.lo, y.hi, yt, s.set.grid, th);
    if (s.result) s.content->paint(fb, p, *s.result, done, y, s.set.thickness, th);
    fb.rect_outline({p.x - 1, p.y - 1, p.w + 2, p.h + 2}, th.fg);

    s.legend_hits.clear();
    if (s.set.legend) s.legend_hits = draw_legend(fb, p, s.content->legend(th), th);
    draw_title(fb, s.layout, s.set.title, th);
    if (mode_label) fb.text(Framebuffer::kCharW, s.layout.readout.y + 4, cmode_name(s.set.cmode), th.fg);
}

void render_headless(Framebuffer& fb, Screen& s, const Theme& th) {
    s.layout = compute_layout(fb.width(), fb.height());
    s.result = s.content->new_result(s.xview(), s.set.index, s.views.top().y);
    s.content->reduce(*s.result, {});
    s.result->ended = true;
    compose(fb, s, s.result->width(), th, /*mode_label=*/true);
}

}  // namespace ssp
