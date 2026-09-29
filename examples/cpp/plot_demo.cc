// Interactive xplot demo: three complex traces with different xstart/xdelta (like the old
// tests/test_plot.cc), sized so the first render is visibly progressive.
//   bazel run -c opt //examples/cpp:plot_demo -- [samples_per_trace]
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "ssp/plot.h"

int main(int argc, char** argv) {
    const std::size_t n = argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 50'000'000;
    using C = std::complex<float>;
    std::vector<C> a(n), b(n), c(n / 4);
    for (std::size_t i = 0; i < n; ++i) {
        const float t = i * 0.02f;
        a[i] = C(std::sin(2 * float(M_PI) * t / 1e3f), 0.0f);
        b[i] = C(0.0f, std::sin(2 * float(M_PI) * i * 0.01f / 1e3f) + 1.0f);
    }
    for (std::size_t i = 0; i < c.size(); ++i) c[i] = std::polar(1.0f, 2 * float(M_PI) * i * 0.02f / 1e3f);

    ssp::Signal sb(b, 110.0, 0.01, "Imag Sin + 1");
    sb.style = ssp::Style::Dots;
    ssp::PlotOptions o;
    o.title = "ssp demo - " + std::to_string(n) + " samples/trace";
    o.thickness = 1;
    try {
        ssp::plot({ssp::Signal(a, 100.0, 0.02, "Real Sin"), sb,
                   ssp::Signal(c, 100.0 + 0.75 * n * 0.02, 0.02, "CPLX")},
                  o);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "ssp_demo: %s\n", e.what());
        return 1;
    }
}
