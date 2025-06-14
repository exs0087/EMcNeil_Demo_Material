#include "ode4.h"

std::vector<Vec7> ode4(
    const std::function<Vec7(double, const Vec7&)>& f,
    const std::vector<double>& t,
    const Vec7& y0
) {
    size_t N = t.size();
    std::vector<Vec7> Y(N);
    Y[0] = y0;

    for (size_t i = 1; i < N; ++i) {
        double h = t[i] - t[i-1];
        const Vec7& y_prev = Y[i-1];
        double t_prev = t[i-1];

        auto k1 = f(t_prev, y_prev);
        auto k2 = f(t_prev + 0.5*h, y_prev + k1*(0.5*h));
        auto k3 = f(t_prev + 0.5*h, y_prev + k2*(0.5*h));
        auto k4 = f(t[i],       y_prev + k3*h);

        Y[i] = y_prev + (k1 + k2*2.0 + k3*2.0 + k4)*(h/6.0);
    }

    return Y;
}
