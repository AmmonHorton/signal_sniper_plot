/// @file raster_content.h
/// @brief xraster data: frames of samples shown as an image, plus the row/column views that
/// x- and y-cuts plot.
#pragma once

#include <optional>

#include "app/content.h"
#include "ssp/plot.h"

namespace ssp {

class RasterContent final : public Content {
public:
    /// @throws std::invalid_argument for unusable data or options.
    RasterContent(Signal data, RasterOptions options);

    bool is_raster() const override { return true; }
    Settings initial_settings() const override;
    View home(const Settings& s) const override;
    bool supports(CMode m) const override { return m != CMode::Auto && m != CMode::IR; }
    Range xaxis() const override { return {sig_.xstart, sig_.xdelta}; }

    std::shared_ptr<ReduceResult> new_result(const ReduceRequest& req) const override;
    bool reusable(const ReduceResult& r, const ReduceRequest& req) const override;
    void reduce(ReduceResult& r, const CancelToken& ct,
                const std::function<void()>& progress = {}) const override;
    /// @brief z extent of rows [0, done).
    Span extent(const ReduceResult& r, int done) const override;
    void paint(Framebuffer& fb, const ReduceResult& r, const PaintArgs& a) const override;
    std::optional<double> z_at(double x, double y, const Settings& s) const override;

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return opts_.subsize; }
    const RasterOptions& options() const { return opts_; }

    // The only place cells are addressed (so a ring-buffer waterfall only has to change these).
    /// @brief Cells [c0, c1) of row r as a trace with the column x axis.
    Signal row_signal(std::size_t r, std::size_t c0, std::size_t c1) const;
    /// @brief Cells [r0, r1) of column c as a trace with the frame (y) axis.
    Signal col_signal(std::size_t c, std::size_t r0, std::size_t r1) const;
    /// @brief Cells of row r that exist (the last frame may be partial).
    std::size_t cols_in_row(std::size_t r) const;

    /// @brief Column/row holding data coordinate x/y, if any.
    std::optional<std::size_t> col_at(double x, bool index) const;
    std::optional<std::size_t> row_at(double y, bool index) const;
    /// @brief Columns/rows whose cells overlap [lo, hi] (a zoom box), as [first, last).
    SampleRange cols_in(Range x, bool index) const;
    SampleRange rows_in(Range y, bool index) const;

    /// @brief {start, delta} of the column (x) and row (y) axes.
    Range col_axis(bool index) const { return index ? Range{0.0, 1.0} : Range{sig_.xstart, sig_.xdelta}; }
    Range row_axis(bool index) const { return index ? Range{0.0, 1.0} : Range{opts_.ystart, opts_.ydelta}; }

private:
    Signal sig_;
    RasterOptions opts_;
    std::size_t rows_ = 0;
    CMode initial_cmode_ = CMode::Real;
};

}  // namespace ssp
