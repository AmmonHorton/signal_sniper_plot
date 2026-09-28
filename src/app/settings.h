/// @file settings.h
/// @brief Per-screen user settings and the zoom stack.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "ssp/types.h"

namespace ssp {

struct ReduceResult;

/// @brief What the x readout shows: the x value, the sample index, or 1/x (SigPlot 'A').
enum class Absc : uint8_t { X, Index, Inverse };

/// @brief What the user can change at run time (keys, legend). Owned by one screen.
struct Settings {
    std::string title;
    CMode cmode = CMode::Real;  ///< Always resolved (never Auto).
    PhaseUnits phunits = PhaseUnits::Radians;
    bool index = false;
    bool grid = true;
    bool legend = true;
    bool cross = true;  ///< Crosshair over the data area.
    Absc absc = Absc::X;
    // Raster only.
    Colormap cmap = Colormap::Ramp;
    Reduce reduce = Reduce::Max;
    std::optional<Range> zfixed;  ///< Fixed colour range; empty = autoscale.
    int thickness = 1;
};

/// @brief One zoom level. `y` empty = autoscale y from the data in view.
struct View {
    Range x;
    std::optional<Range> y;
    /// Made by the wheel or panning: further wheel/pan steps change this level in place
    /// instead of stacking up new ones.
    bool adjustable = false;
    /// IR home: x (the real part) autoscales from the data, so `x` is not used.
    bool auto_x = false;
    /// Last reduce run for this level, so unzooming back to it can reuse it instead of
    /// recomputing. Reused only if complete and still matching (TraceContent::reusable).
    std::shared_ptr<ReduceResult> cache;
};

/// @brief SigPlot's Mx.stk: level 0 is "home", zoom pushes, unzoom pops.
class ViewStack {
public:
    static constexpr std::size_t kMaxLevels = 10;

    void reset(View home) { lv_.assign(1, home); }
    const View& top() const { return lv_.back(); }
    std::size_t level() const { return lv_.size() - 1; }  ///< 0 = home

    /// @return false (and nothing changes) when already kMaxLevels deep.
    bool push(View v) {
        if (lv_.size() >= kMaxLevels) return false;
        lv_.push_back(v);
        return true;
    }
    /// @return false when already at home.
    bool pop() {
        if (lv_.size() <= 1) return false;
        lv_.pop_back();
        return true;
    }
    void pop_all() { lv_.resize(1); }
    View& home() { return lv_.front(); }
    /// @brief Make every level autoscale y (used when the plotted quantity changes).
    void clear_y() {
        for (auto& v : lv_) v.y.reset();
    }
    /// @brief Forget every level's cached result (their inputs changed).
    void clear_cache() {
        for (auto& v : lv_) v.cache.reset();
    }
    View& top_mut() { return lv_.back(); }
    /// @brief Level below the top (the one right-click returns to); home when at home.
    const View& parent() const { return lv_.size() > 1 ? lv_[lv_.size() - 2] : lv_.front(); }

private:
    std::vector<View> lv_{View{{-1.0, 1.0}, std::nullopt}};
};

}  // namespace ssp
