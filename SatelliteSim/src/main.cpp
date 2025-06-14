// src/main.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <filesystem>
#include <cmath>

#include "adcs_hal.h"       // HAL_init, HAL_set_control_mode(), ControlMode
#include "Vec3.h"
#include "Vec4.h"
#include "Vec7.h"
#include "eulerseqns2.h"    // eulerseqns2(), getActHistory()
#include "ode4.h"
#include "magnetic_field.h" // getBodyFieldAt()

// Where to drop our CSV (set via CMake with -DOUTPUT_DIR="..."):
#ifndef OUTPUT_DIR
#error "OUTPUT_DIR must be defined by CMake"
#endif
constexpr char const* kOutputDir = OUTPUT_DIR;

/// One row in the CSV
struct Row {
    double time;
    double wx, wy, wz;
    double q1, q2, q3, q4;
    double Bx, By, Bz;
    double mx, my, mz;
};

int main(int argc, char** argv) {
    //
    // 1) Simulation parameters
    //
    const double sampleRate = 1.0;    // seconds
    const double t_final    = 16000.0; // seconds
    const size_t N = static_cast<size_t>(t_final / sampleRate) + 1;

    //
    // 2) Initialize HAL (tell the embedded stub our control‐step interval)
    //
    HAL_init(sampleRate);

    // Optional: let user pick the starting mode on the command line
    if (argc > 1 && std::string(argv[1]) == "bdot") {
        HAL_set_control_mode(CTRL_MODE_B_DOT);
        std::cout << "Starting in Ḃ‐law (bang–bang) mode\n";
    } else {
        HAL_set_control_mode(CTRL_MODE_WX_B);
        std::cout << "Starting in ω×B proportional mode\n";
    }

    //
    // 3) Build time array
    //
    std::vector<double> times(N);
    for (size_t i = 0; i < N; ++i) {
        times[i] = i * sampleRate;
    }

    //
    // 4) Initial state [ωx,ωy,ωz, qx,qy,qz,qw]
    //
    Vec7 y0{ 0.0, 0.0, 0.0,   0.0, 0.0, 0.0, 1.0 };

    //
    // 5) Integrate dynamics (RK4), timing the run
    //
    auto t_start = std::chrono::high_resolution_clock::now();
    auto Y       = ode4(eulerseqns2, times, y0);
    auto t_end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << "Simulation took " << elapsed.count() << " seconds.\n";

    //
    // 6) Grab the one‐per‐step actuator commands we recorded in eulerseqns2()
    //
    const auto& act_hist = getActHistory();

    //
    // 7) Assemble one Row per step
    //
    std::vector<Row> allResults;
    allResults.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        Row r;
        r.time = times[i];

        // body‐rates
        r.wx = Y[i][0];
        r.wy = Y[i][1];
        r.wz = Y[i][2];

        // quaternion
        r.q1 = Y[i][3];
        r.q2 = Y[i][4];
        r.q3 = Y[i][5];
        r.q4 = Y[i][6];

        // recompute body‐frame B via **seconds-based** sim HAL
        Vec4 qt{ r.q1, r.q2, r.q3, r.q4 };
        Vec3 B = getBodyFieldAt(times[i], qt);
        r.Bx = B.x;
        r.By = B.y;
        r.Bz = B.z;

        // actual dipole commands from HAL
        r.mx = act_hist[i].mx;
        r.my = act_hist[i].my;
        r.mz = act_hist[i].mz;

        allResults.push_back(r);
    }

    //
    // 8) Write out CSV for analysis
    //
    namespace fs = std::filesystem;
    fs::create_directories(kOutputDir);  // ensure the folder exists

    std::ofstream csv{ fs::path(kOutputDir) / "sim_results.csv" };
    csv << "time,wx,wy,wz,"
           "q1,q2,q3,q4,"
           "Bx,By,Bz,"
           "mx,my,mz\n";

    // Helper: print a number or blank if non-finite
    auto writeOrBlank = [&](double x) {
        if (std::isfinite(x)) csv << x;
        // else leave empty
    };

    for (auto const& r : allResults) {
        writeOrBlank(r.time); csv << ',';
        writeOrBlank(r.wx);   csv << ',';
        writeOrBlank(r.wy);   csv << ',';
        writeOrBlank(r.wz);   csv << ',';
        writeOrBlank(r.q1);   csv << ',';
        writeOrBlank(r.q2);   csv << ',';
        writeOrBlank(r.q3);   csv << ',';
        writeOrBlank(r.q4);   csv << ',';
        writeOrBlank(r.Bx);   csv << ',';
        writeOrBlank(r.By);   csv << ',';
        writeOrBlank(r.Bz);   csv << ',';
        writeOrBlank(r.mx);   csv << ',';
        writeOrBlank(r.my);   csv << ',';
        writeOrBlank(r.mz);
        csv << '\n';
    }
    csv.close();

    std::cout << "Wrote sim results to "
              << (fs::path(kOutputDir) / "sim_results.csv") << "\n";

    return 0;
}
