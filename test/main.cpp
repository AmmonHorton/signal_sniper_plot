#include "plot.hpp"
#include <vector>
#include <cmath>

int main() {
    constexpr size_t N = 1e4;
    std::vector<float> signal;
    for (size_t i = 0; i < N; ++i) {
        float t = i * 0.01f;
        signal.push_back(std::sin(2 * M_PI * t));
    }

    xplot::plot_buffer(
        signal.data(),
        signal.size(),
        sizeof(float),
        false,  // not complex
        true,   // is float
        2524608000.0,  // J1950 start time (2020)
        0.01,           // xdelta (s)
        "Signal Plot - Frame 37",
        2,
        std::make_optional<std::pair<double, double>>({-5, 5})
    );
    
    return 0;
}