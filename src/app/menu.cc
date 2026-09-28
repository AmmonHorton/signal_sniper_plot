#include "app/menu.h"

#include <algorithm>

#include "app/controller.h"

namespace ssp {
namespace {

constexpr int kCW = Framebuffer::kCharW;
constexpr int kItemH = Framebuffer::kCharH + 3;
constexpr int kPad = 3;
constexpr uint32_t kEnter = 0xff0d, kKpEnter = 0xff8d;

MenuItem toggle(const char* label, const char* key, bool on, std::function<unsigned(Screen&)> f) {
    return {label, key, on, std::move(f), {}};
}

struct Box {
    Rect rect;
    const std::vector<MenuItem>* items;
};

int box_width(const std::vector<MenuItem>& items) {
    std::size_t label = 0, key = 0;
    bool arrows = false;
    for (const auto& it : items) {
        label = std::max(label, it.label.size());
        key = std::max(key, it.key.size());
        arrows = arrows || !it.sub.empty();
    }
    const std::size_t chars = 2 + label + (key ? 3 + key : 0) + (arrows ? 2 : 0) + 1;
    return static_cast<int>(chars) * kCW + 2 * kPad;
}

Rect item_rect(const Box& b, int i) {
    return {b.rect.x + 1, b.rect.y + kPad + i * kItemH, b.rect.w - 2, kItemH};
}

/// Boxes of the open levels, each clamped into the window.
std::vector<Box> boxes(const Screen& s, const std::vector<MenuItem>& root) {
    std::vector<Box> out;
    const MenuState& m = *s.ui.menu;
    const std::vector<MenuItem>* items = &root;
    int x = m.x, y = m.y;
    for (std::size_t k = 0; k < m.sel.size(); ++k) {
        Rect r{x, y, box_width(*items), static_cast<int>(items->size()) * kItemH + 2 * kPad};
        r.x = std::max(0, std::min(r.x, s.window_w() - r.w));
        r.y = std::max(0, std::min(r.y, s.window_h() - r.h));
        out.push_back({r, items});
        const int sel = m.sel[k];
        if (k + 1 == m.sel.size() || sel < 0) break;
        x = r.right() - 2;
        y = r.y + sel * kItemH;
        items = &(*items)[sel].sub;
    }
    return out;
}

/// Deepest open level and item under (x, y); level -1 if none.
std::pair<int, int> hit(const std::vector<Box>& bx, int x, int y) {
    for (int k = static_cast<int>(bx.size()) - 1; k >= 0; --k) {
        for (int i = 0; i < static_cast<int>(bx[k].items->size()); ++i) {
            if (item_rect(bx[k], i).contains(x, y)) return {k, i};
        }
        if (bx[k].rect.contains(x, y)) return {k, -1};
    }
    return {-1, -1};
}

const MenuItem* item_at(const std::vector<Box>& bx, int level, int i) {
    if (level < 0 || i < 0 || level >= static_cast<int>(bx.size())) return nullptr;
    return &(*bx[level].items)[i];
}

/// Highlight item i of `level`, opening its submenu if it has one.
void highlight(MenuState& m, const std::vector<Box>& bx, int level, int i) {
    m.sel.resize(level + 1);
    m.sel[level] = i;
    const MenuItem* it = item_at(bx, level, i);
    if (it && !it->sub.empty()) m.sel.push_back(-1);
}

unsigned run(Screen& s, const MenuItem& it) {
    s.ui.menu.reset();
    return (it.action ? it.action(s) : kNone) | kOverlay;
}

}  // namespace

std::vector<MenuItem> build_menu(const Screen& s) {
    auto mode = [&](const char* label, const char* key, CMode m) {
        return toggle(label, key, s.set.cmode == m, [m](Screen& sc) { return action::set_mode(sc, m); });
    };
    auto units = [&](const char* label, PhaseUnits u) {
        return toggle(label, "", s.set.phunits == u, [u](Screen& sc) { return action::set_phunits(sc, u); });
    };
    std::vector<MenuItem> traces;
    for (std::size_t t = 0; t < s.content->size(); ++t) {
        traces.push_back(toggle(s.content->signal(t).name.c_str(), "", s.content->signal(t).visible,
                                [t](Screen& sc) { return action::toggle_trace(sc, t); }));
    }
    return {
        {"Mode", "", false, {},
         {mode("Magnitude", "1", CMode::Mag), mode("Phase", "2", CMode::Phase),
          mode("Real", "3", CMode::Real), mode("Imaginary", "4", CMode::Imag),
          mode("10*log10", "6", CMode::Log10), mode("20*log10", "7", CMode::Log20)}},
        {"Phase units", "", false, {},
         {units("Radians", PhaseUnits::Radians), units("Degrees", PhaseUnits::Degrees),
          units("Cycles", PhaseUnits::Cycles)}},
        {"Scaling", "", false, {},
         {toggle("X range...", "", false,
                 [](Screen& sc) { return action::open_prompt(sc, PromptState::Kind::XRange); }),
          toggle("Y range...", "", false,
                 [](Screen& sc) { return action::open_prompt(sc, PromptState::Kind::YRange); }),
          toggle("Autoscale Y", "", false, action::autoscale_y),
          toggle("Unzoom all", "Home", false, action::unzoom_all)}},
        {"Traces", "", false, {}, traces},
        {},
        toggle("Grid", "g", s.set.grid, [](Screen& sc) { sc.set.grid = !sc.set.grid; return unsigned{kRepaint}; }),
        toggle("Legend", "l", s.set.legend,
               [](Screen& sc) { sc.set.legend = !sc.set.legend; return unsigned{kRepaint}; }),
        toggle("Crosshair", "", s.set.cross, [](Screen& sc) { sc.set.cross = !sc.set.cross; return unsigned{kOverlay}; }),
        toggle("Index x-axis", "i", s.set.index, action::toggle_index),
        {},
        toggle("Save PNG", "Ctrl-S", false, [](Screen&) { return unsigned{kSave}; }),
        toggle("Keypress help", "?", false, [](Screen& sc) { sc.ui.help = true; return unsigned{kOverlay}; }),
        toggle("Exit", "q", false, [](Screen&) { return unsigned{kQuit}; }),
    };
}

unsigned open_menu(Screen& s, int x, int y) {
    s.ui.dragging = false;
    s.ui.menu = MenuState{x + 2, y + 2, {-1}};
    return kOverlay;
}

unsigned handle_menu_event(Screen& s, const InputEvent& e) {
    using T = InputEvent::Type;
    const std::vector<MenuItem> root = build_menu(s);
    const std::vector<Box> bx = boxes(s, root);
    MenuState& m = *s.ui.menu;

    switch (e.type) {
        case T::Motion: {
            const auto [level, i] = hit(bx, e.x, e.y);
            const MenuItem* it = item_at(bx, level, i);
            if (it && !it->separator()) highlight(m, bx, level, i);
            return kOverlay;
        }
        case T::Press:
        case T::Release: {
            const auto [level, i] = hit(bx, e.x, e.y);
            const MenuItem* it = item_at(bx, level, i);
            if (e.type == T::Press && level < 0) {  // click outside closes, and is consumed
                s.ui.menu.reset();
                return kOverlay;
            }
            if (!it || it->separator()) return kNone;
            if (!it->sub.empty()) {
                if (e.type == T::Press) highlight(m, bx, level, i);
                return kOverlay;
            }
            // Click an action, or middle-button press-drag-release onto one.
            if (e.type == T::Press ? e.button != 2 : e.button == 2) return run(s, *it);
            return kNone;
        }
        case T::Key: {
            const int level = static_cast<int>(m.sel.size()) - 1;
            const auto& items = *bx[std::min<std::size_t>(level, bx.size() - 1)].items;
            const int n = static_cast<int>(items.size());
            int& cur = m.sel[level];
            if (e.key == key::kEscape) {
                s.ui.menu.reset();
            } else if (e.key == key::kDown || e.key == key::kUp) {
                const int step = e.key == key::kDown ? 1 : -1;
                int next = cur < 0 ? (step > 0 ? -1 : n) : cur;
                do {
                    next = (next + step + n) % n;
                } while (items[next].separator());
                cur = next;
            } else if (e.key == key::kLeft) {
                if (m.sel.size() > 1) m.sel.pop_back();
            } else if (cur >= 0 && (e.key == key::kRight || e.key == kEnter || e.key == kKpEnter ||
                                    e.key == key::kSpace)) {
                const MenuItem& it = items[cur];
                if (!it.sub.empty()) {
                    m.sel.push_back(0);
                } else if (e.key != key::kRight) {
                    return run(s, it);
                }
            }
            return kOverlay;
        }
        default:
            return kNone;
    }
}

std::vector<Rect> draw_menu(Framebuffer& fb, const Screen& s, const Theme& th) {
    std::vector<Rect> touched;
    if (!s.ui.menu) return touched;
    const std::vector<MenuItem> root = build_menu(s);
    const std::vector<Box> bx = boxes(s, root);
    for (std::size_t k = 0; k < bx.size(); ++k) {
        const Box& b = bx[k];
        fb.fill_rect(b.rect, th.bg);
        fb.rect_outline(b.rect, th.fg);
        const int sel = s.ui.menu->sel[k];
        for (int i = 0; i < static_cast<int>(b.items->size()); ++i) {
            const MenuItem& it = (*b.items)[i];
            const Rect r = item_rect(b, i);
            if (it.separator()) {
                fb.hline(r.x + kPad, r.right() - kPad - 1, r.y + kItemH / 2, th.dim);
                continue;
            }
            const bool hot = i == sel;
            const uint32_t fg = hot ? th.bg : th.fg;
            if (hot) fb.fill_rect(r, th.fg);
            const int ty = r.y + 2;
            if (it.checked) fb.text(r.x + kPad, ty, "*", fg);
            fb.text(r.x + kPad + 2 * kCW, ty, it.label, fg);
            int right = r.right() - kPad;
            if (!it.sub.empty()) {
                right -= kCW;
                fb.text(right, ty, ">", fg);
                right -= kCW;
            }
            if (!it.key.empty()) {
                fb.text(right - Framebuffer::text_width(it.key), ty, it.key, hot ? th.bg : th.dim);
            }
        }
        touched.push_back(b.rect);
    }
    return touched;
}

}  // namespace ssp
