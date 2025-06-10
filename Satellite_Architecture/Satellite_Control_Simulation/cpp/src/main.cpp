// main.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include "ode4m.h"
#include "eulerseqns2.h"

// Declare the global MHIST from ode4m.cpp so we can write it out here.
// In ode4m.cpp, MHIST must not be static for this to link correctly.
extern std::vector<std::array<double,3>> MHIST;

int main() {
    // 1) Build tspan = [0, 1, 2, ..., 16000]
    std::vector<double> tspan(16001);
    for (int i = 0; i <= 16000; ++i) {
        tspan[i] = static_cast<double>(i);
    }

    // 2) Initial state: w0 = [0, 0.2, 0.2], q0 = [0, 0, 0, 1]
    //    State vector y0 = { w1, w2, w3, q1, q2, q3, q4 }
    std::array<double,7> y0 = { 0.0, 0.2, 0.2,   0.0, 0.0, 0.0, 1.0 };

    // Define a lambda that simply calls eulerseqns2
    auto solver = [](double t, const std::array<double,7>& y) {
        return eulerseqns2(t, y);
    };

    // -------------------
    // 3) Run for cmode = 1 (B-dot)
    // -------------------
    int cmode1 = 1;
    auto Ysim1 = ode4m(solver, tspan, y0, cmode1);

    // 4) Save Ysim1 to "ysim_cmode1.csv"
    {
        std::ofstream ofs("ysim_cmode1.csv");
        ofs << "time,w1,w2,w3,q1,q2,q3,q4\n";
        for (size_t i = 0; i < Ysim1.size(); ++i) {
            double t = tspan[i];
            const auto &y = Ysim1[i];
            ofs << t << ","
                << y[0] << "," << y[1] << "," << y[2] << ","
                << y[3] << "," << y[4] << "," << y[5] << "," << y[6] << "\n";
        }
        ofs.close();
    }

    // 5) Save MHIST for cmode = 1 to "mhist_cmode1.csv"
    {
        std::ofstream ofs("mhist_cmode1.csv");
        ofs << "time,m1,m2,m3\n";
        for (size_t i = 0; i < MHIST.size(); ++i) {
            double t = tspan[i];
            const auto &m = MHIST[i];
            ofs << t << "," << m[0] << "," << m[1] << "," << m[2] << "\n";
        }
        ofs.close();
    }

    // 6) Compute and print final value metric for cmode = 1
    {
        // Last state → extract ω components
        auto y_end = Ysim1.back();
        double wf1 = y_end[0], wf2 = y_end[1], wf3 = y_end[2];
        double wf_mag = std::sqrt(wf1*wf1 + wf2*wf2 + wf3*wf3);
        double mu = 398600.0;    // km^3/s^2
        double Re = 6371.0;      // km
        double altitude = 400.0; // km
        double r = Re + altitude;
        double Torbit = 2.0*M_PI * std::sqrt((r*r*r)/mu);
        double final_metric1 = wf_mag * (Torbit / (2.0*M_PI));
        std::cout << "B-dot final value (rev/orbit): " << final_metric1 << std::endl;
    }

    // -------------------
    // 7) Run for cmode = 2 (ω×B)
    // -------------------
    int cmode2 = 2;
    auto Ysim2 = ode4m(solver, tspan, y0, cmode2);

    // 8) Save Ysim2 to "ysim_cmode2.csv"
    {
        std::ofstream ofs("ysim_cmode2.csv");
        ofs << "time,w1,w2,w3,q1,q2,q3,q4\n";
        for (size_t i = 0; i < Ysim2.size(); ++i) {
            double t = tspan[i];
            const auto &y = Ysim2[i];
            ofs << t << ","
                << y[0] << "," << y[1] << "," << y[2] << ","
                << y[3] << "," << y[4] << "," << y[5] << "," << y[6] << "\n";
        }
        ofs.close();
    }

    // 9) Save MHIST for cmode = 2 to "mhist_cmode2.csv"
    {
        std::ofstream ofs("mhist_cmode2.csv");
        ofs << "time,m1,m2,m3\n";
        for (size_t i = 0; i < MHIST.size(); ++i) {
            double t = tspan[i];
            const auto &m = MHIST[i];
            ofs << t << "," << m[0] << "," << m[1] << "," << m[2] << "\n";
        }
        ofs.close();
    }

    // 10) Compute and print final value metric for cmode = 2
    {
        auto y_end = Ysim2.back();
        double wf1 = y_end[0], wf2 = y_end[1], wf3 = y_end[2];
        double wf_mag2 = std::sqrt(wf1*wf1 + wf2*wf2 + wf3*wf3);
        double mu = 398600.0, Re = 6371.0, altitude = 400.0;
        double r = Re + altitude;
        double Torbit = 2.0*M_PI * std::sqrt((r*r*r)/mu);
        double final_metric2 = wf_mag2 * (Torbit / (2.0*M_PI));
        std::cout << "ω×B final value (rev/orbit): " << final_metric2 << std::endl;
    }

    return 0;
}
