#include "ssp/app/screen.h"

#include "ssp/app/raster_content.h"

namespace ssp {

Screen::Screen(std::unique_ptr<Content> c) : content(std::move(c)) {
    set = content->initial_settings();
    if (TraceContent* t = traces(); t && set.cmode == CMode::IR) ir = t->initial_ir();
    go_home();
}

Range Screen::x_for(int done) const {
    if (!views.top().auto_x) return views.top().x;
    return (result && done > 0 && result->x) ? *result->x : xshown;
}

const RasterContent* Screen::raster() const {
    return dynamic_cast<const RasterContent*>(content.get());
}

Range Screen::z_for(int done) const {
    if (set.zfixed) return *set.zfixed;
    if (!result) return zshown;
    const Span e = content->extent(*result, done);
    if (e.empty()) return zshown;
    return e.hi > e.lo ? Range{e.lo, e.hi} : Range{e.lo - 1.0, e.hi + 1.0};
}

Range Screen::level_x(const View& v) const {
    if (!v.auto_x) return v.x;
    return (v.cache && v.cache->complete() && v.cache->x) ? *v.cache->x : xshown;
}

ReduceRequest Screen::request() const {
    const View& v = views.top();
    ReduceRequest q;
    q.view = xview();
    q.index = set.index;
    q.auto_x = v.auto_x;
    q.y = v.y;
    q.reduce = set.reduce;
    if (set.cmode == CMode::IR) q.ir = ir;
    return q;
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
    const Range x = xshown;
    return x.lo + (px - layout.plot.x + 0.5) / layout.plot.w * (x.hi - x.lo);
}

double Screen::py_to_y(int py) const {
    const double f = (py - layout.plot.y + 0.5) / layout.plot.h * (yshown.hi - yshown.lo);
    return ydown() ? yshown.lo + f : yshown.hi - f;
}

}  // namespace ssp
