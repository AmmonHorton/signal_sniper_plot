#include "ssp/render/framebuffer.h"

#include <gtest/gtest.h>

#include <cmath>

#include "ssp/render/png.h"

namespace ssp {
namespace {

int count_lit(const Framebuffer& fb) {
    int n = 0;
    for (int y = 0; y < fb.height(); ++y) {
        for (int x = 0; x < fb.width(); ++x) n += fb.at(x, y) != 0;
    }
    return n;
}

TEST(Framebuffer, LineIsClippedToRect) {
    Framebuffer fb(20, 20);
    const Rect clip{5, 5, 10, 10};
    fb.line(-100, 10.5, 100, 10.5, 1, 1, clip);
    EXPECT_EQ(count_lit(fb), 10);
    for (int x = 5; x < 15; ++x) EXPECT_EQ(fb.at(x, 10), 1u);
}

TEST(Framebuffer, LineEntirelyOutsideDrawsNothing) {
    Framebuffer fb(20, 20);
    fb.line(-5, -5, -1, 30, 1, 3, fb.bounds());
    fb.line(0, 0, 1e300, std::nan(""), 1, 1, fb.bounds());
    EXPECT_EQ(count_lit(fb), 0);
}

TEST(Framebuffer, ThickLineUsesSquareBrush) {
    Framebuffer fb(20, 20);
    fb.line(2.5, 10.5, 17.5, 10.5, 1, 3, fb.bounds());
    EXPECT_EQ(count_lit(fb), 18 * 3);  // columns 2..17 plus the brush overhang at each end
}

TEST(Framebuffer, TextDrawsGlyphs) {
    Framebuffer fb(40, 20);
    fb.text(1, 1, "A", 7);
    EXPECT_GT(count_lit(fb), 10);
    EXPECT_EQ(Framebuffer::text_width("abc"), 18);
}

TEST(Png, RoundTrips) {
    Framebuffer fb(17, 9, 0x123456);
    fb.line(0, 0, 16, 8, 0xFF8000, 1, fb.bounds());
    const Framebuffer back = decode_png(encode_png(fb));
    ASSERT_EQ(back.width(), 17);
    ASSERT_EQ(back.height(), 9);
    for (int y = 0; y < 9; ++y) {
        for (int x = 0; x < 17; ++x) ASSERT_EQ(back.at(x, y), fb.at(x, y));
    }
}

}  // namespace
}  // namespace ssp
