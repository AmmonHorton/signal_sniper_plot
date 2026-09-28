#include "render/axis.h"

#include <gtest/gtest.h>

namespace ssp {
namespace {

using Labels = std::vector<std::string>;

TEST(NiceTics, StepsAreOneTwoTwoPointFiveFive) {
    EXPECT_DOUBLE_EQ(nice_tics(0, 10, 5).dtic, 2.0);
    EXPECT_DOUBLE_EQ(nice_tics(0, 10, 4).dtic, 2.5);
    EXPECT_DOUBLE_EQ(nice_tics(0, 10, 2).dtic, 5.0);
    EXPECT_DOUBLE_EQ(nice_tics(0, 1e-6, 5).dtic, 2e-7);
    EXPECT_DOUBLE_EQ(nice_tics(-3.3, 7.1, 5).dtic, 2.0);
}

TEST(MakeTicks, PlainRange) {
    const AxisTicks t = make_ticks(0, 1, 5);
    EXPECT_EQ(t.labels, (Labels{"0.0", "0.2", "0.4", "0.6", "0.8", "1.0"}));
    EXPECT_EQ(t.note, "");
}

TEST(MakeTicks, TwoPointFiveStepGetsADecimal) {
    EXPECT_EQ(make_ticks(0, 10, 4).labels, (Labels{"0.0", "2.5", "5.0", "7.5", "10.0"}));
}

TEST(MakeTicks, NegativeZeroIsPrintedAsZero) {
    const AxisTicks t = make_ticks(-1, 1, 4);
    EXPECT_EQ(t.labels, (Labels{"-1.0", "-0.5", "0.0", "0.5", "1.0"}));
}

TEST(MakeTicks, TinyXdeltaUsesEngineeringMultiplier) {
    // Regression: "%.2f" printed every tick as 0.00.
    const AxisTicks t = make_ticks(0, 1e-6, 5);
    EXPECT_EQ(t.note, "x1e-6");
    EXPECT_EQ(t.labels, (Labels{"0.0", "0.2", "0.4", "0.6", "0.8", "1.0"}));
}

TEST(MakeTicks, LargeValuesUseEngineeringMultiplier) {
    const AxisTicks t = make_ticks(0, 5e4, 5);
    EXPECT_EQ(t.note, "x1e3");
    EXPECT_EQ(t.labels.back(), "50");
}

TEST(MakeTicks, HugeOffsetSwitchesToRelativeLabels) {
    const AxisTicks t = make_ticks(1e9, 1e9 + 10, 5);
    EXPECT_EQ(t.note, "+1000000000");
    EXPECT_EQ(t.labels, (Labels{"0", "2", "4", "6", "8", "10"}));
}

TEST(MakeTicks, DegenerateRangeHasNoTicks) {
    EXPECT_TRUE(make_ticks(1, 1, 5).values.empty());
}

}  // namespace
}  // namespace ssp
