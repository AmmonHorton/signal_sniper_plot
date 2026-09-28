#include "app/cut.h"

#include <algorithm>

#include "app/raster_content.h"
#include "render/axis.h"

namespace ssp {
namespace {

struct Built {
    Signal signal;
    std::string title;
};

Built build(const CutInfo& c) {
    const std::string& raster_title = c.raster_title;
    const RasterContent& rc = *c.source;
    Built b;
    std::string where;
    if (c.xcut) {
        const std::size_t last = std::min(c.span.second, rc.cols_in_row(c.index));
        b.signal = rc.row_signal(c.index, c.span.first, std::max(last, c.span.first));
        const Range ry = rc.row_axis(c.index_axes);
        where = "row " + std::to_string(c.index) + " (y = " + format_g(ry.lo + c.index * ry.hi, 6) + ")";
    } else {
        std::size_t last = c.span.second;  // a partial last frame may not reach this column
        while (last > c.span.first && rc.cols_in_row(last - 1) <= c.index) --last;
        b.signal = rc.col_signal(c.index, c.span.first, last);
        const Range cx = rc.col_axis(c.index_axes);
        where = "column " + std::to_string(c.index) + " (x = " + format_g(cx.lo + c.index * cx.hi, 6) + ")";
    }
    if (c.index_axes) {  // cells labelled by number, like the raster's axis
        b.signal.xstart = static_cast<double>(c.span.first);
        b.signal.xdelta = 1.0;
    }
    b.signal.name = (c.xcut ? "row " : "column ") + std::to_string(c.index);
    b.title = (raster_title.empty() ? std::string() : raster_title + " - ") +
              (c.xcut ? "x-cut " : "y-cut ") + where;
    return b;
}

}  // namespace

std::unique_ptr<Screen> make_cut(const Screen& raster, bool xcut, std::size_t index) {
    const RasterContent& rc = *raster.raster();
    const bool idx = raster.set.index;
    const SampleRange cols = rc.cols_in(raster.xshown, idx), rows = rc.rows_in(raster.yshown, idx);
    CutInfo c{&rc, xcut, index, xcut ? rows : cols, xcut ? cols : rows, idx, raster.set.title};
    const Built b = build(c);

    PlotOptions po;
    po.title = b.title;
    po.cmode = raster.set.cmode;
    po.phunits = raster.set.phunits;
    po.legend = false;
    auto s = std::make_unique<Screen>(std::make_unique<TraceContent>(std::vector<Signal>{b.signal}, po));
    s->cut = c;
    s->message = std::string("Esc: back to raster   PgUp/PgDn: previous/next ") + (xcut ? "row" : "column");
    return s;
}

bool step_cut(Screen& cut, int step) {
    CutInfo& c = *cut.cut;
    const long long next = static_cast<long long>(c.index) + step;
    if (next < static_cast<long long>(c.steps.first) || next >= static_cast<long long>(c.steps.second)) {
        cut.message = std::string("no more ") + (c.xcut ? "rows" : "columns") + " in the zoom box";
        return false;
    }
    c.index = static_cast<std::size_t>(next);
    const Built b = build(c);

    PlotOptions po;
    po.cmode = cut.set.cmode;
    po.legend = false;
    cut.content = std::make_unique<TraceContent>(std::vector<Signal>{b.signal}, po);
    cut.set.title = b.title;
    cut.views.clear_cache();  // same zoom, different data
    cut.result.reset();
    return true;
}

}  // namespace ssp
