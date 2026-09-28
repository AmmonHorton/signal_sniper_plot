#include "app/screen.h"

namespace ssp {

Screen::Screen(std::unique_ptr<TraceContent> c, const PlotOptions& o)
    : content(std::move(c)), opts(o) {
    set.title = o.title;
    set.cmode = content->initial_cmode();
    set.phunits = o.phunits;
    set.index = o.index;
    set.grid = o.grid;
    set.legend = o.legend;
    set.thickness = o.thickness;
    go_home();
}

View Screen::home() const {
    View v;
    const bool original_axes = set.index == opts.index;
    v.x = (original_axes && opts.xrange) ? *opts.xrange : content->x_extent(set.index);
    if (opts.yrange && set.cmode == content->initial_cmode()) v.y = *opts.yrange;
    return v;
}

Range Screen::y_for(int done) const {
    if (views.top().y) return *views.top().y;
    if (!result) return yshown;
    if (done > 0 && result->y) return *result->y;  // fixed up front (dots mode)
    const Span e = content->extent(*result, done);
    return e.empty() ? yshown : autoscale_y(e);
}

XView Screen::xview() const {
    const Range x = views.top().x;
    return {x.lo, x.hi, layout.plot.w, layout.plot.h, set.cmode, set.phunits};
}

double Screen::px_to_x(int px) const {
    const Range x = views.top().x;
    return x.lo + (px - layout.plot.x + 0.5) / layout.plot.w * (x.hi - x.lo);
}

double Screen::py_to_y(int py) const {
    return yshown.hi - (py - layout.plot.y + 0.5) / layout.plot.h * (yshown.hi - yshown.lo);
}

}  // namespace ssp
