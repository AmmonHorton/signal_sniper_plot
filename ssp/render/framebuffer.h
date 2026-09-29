/// @file framebuffer.h
/// @brief Client-side 0x00RRGGBB image plus the handful of primitives the plots need.
#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

namespace ssp {

struct Rect {
    int x = 0, y = 0, w = 0, h = 0;

    int right() const { return x + w; }   ///< One past the last column.
    int bottom() const { return y + h; }  ///< One past the last row.
    bool empty() const { return w <= 0 || h <= 0; }
    bool contains(int px, int py) const { return px >= x && px < right() && py >= y && py < bottom(); }
    Rect intersect(const Rect& o) const;
};

class Framebuffer {
public:
    static constexpr int kCharW = 6;   ///< Embedded misc-fixed 6x13 font.
    static constexpr int kCharH = 13;

    Framebuffer() = default;
    Framebuffer(int w, int h, uint32_t fill = 0);

    int width() const { return w_; }
    int height() const { return h_; }
    Rect bounds() const { return {0, 0, w_, h_}; }
    const uint32_t* row(int y) const { return px_.data() + static_cast<std::size_t>(y) * w_; }
    uint32_t* row(int y) { return px_.data() + static_cast<std::size_t>(y) * w_; }
    uint32_t at(int x, int y) const { return row(y)[x]; }

    void fill(uint32_t c) { fill_rect(bounds(), c); }
    void fill_rect(Rect r, uint32_t c);
    void rect_outline(Rect r, uint32_t c);

    /// @brief Inclusive horizontal/vertical runs. `dot` > 0 draws every dot-th pixel only.
    void hline(int x0, int x1, int y, uint32_t c, int dot = 0);
    void vline(int x, int y0, int y1, uint32_t c, int dot = 0);

    /// @brief Line between real-valued pixel coordinates, clipped to `clip`. A point p maps
    /// to pixel floor(p); `thick` is the side of a square brush.
    void line(double x0, double y0, double x1, double y1, uint32_t c, int thick, Rect clip);

    /// @brief Filled disk of radius r (r = 0 is a single pixel), clipped to `clip`.
    void disk(int cx, int cy, int r, uint32_t c, Rect clip);

    /// @brief Text in the embedded font; (x, y) is the top-left of the first cell.
    void text(int x, int y, std::string_view s, uint32_t c);
    static int text_width(std::string_view s) { return kCharW * static_cast<int>(s.size()); }

private:
    void brush(int x, int y, uint32_t c, int thick, const Rect& clip);

    int w_ = 0, h_ = 0;
    std::vector<uint32_t> px_;
};

}  // namespace ssp
