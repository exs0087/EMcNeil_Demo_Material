// src/main.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <filesystem>

#include "adcs_hal.h"         // HAL_init, HAL_set_control_mode, ControlMode
#include "Vec3.h"
#include "Vec4.h"
#include "Vec7.h"
#include "eulerseqns2.h"      // eulerseqns2(), getActHistory()
#include "ode4.h"
#include "magnetic_field.h"   // getBodyFieldAt()

// Where to drop our CSV
constexpr char const* OUTPUT_DIR = "../Test_Setup/cpp";

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
    const double t_final   = 16000.0; // seconds
    const size_t N = static_cast<size_t>(t_final / sampleRate) + 1;

    //
    // 2) Initialize HAL (tell the embedded stub our control‐step interval)
    //
    HAL_init(sampleRate);

    // Optional: let user pick the starting mode on the command line:
    //   argv[1] == "bdot" -> start in Ḃ‐law, else ω×B
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
    auto Y = ode4(eulerseqns2, times, y0);
    auto t_end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << "Simulation took " << elapsed.count() << " seconds.\n";

    //
    // 6) Grab the one‐per‐step actuator commands we recorded in eulerseqns2()
    //
    auto const& act_hist = getActHistory();

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

        // recompute body‐frame B via sim HAL
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
    fs::create_directories(OUTPUT_DIR);  // ensure the folder exists

    std::ofstream csv{ fs::path(OUTPUT_DIR) / "sim_results.csv" };
    csv << "time,wx,wy,wz,q1,q2,q3,q4,"
           "Bx,By,Bz,mx,my,mz\n";
    for (auto const& r : allResults) {
        csv
          << r.time << ','
          << r.wx   << ',' << r.wy << ',' << r.wz << ','
          << r.q1   << ',' << r.q2 << ',' << r.q3 << ',' << r.q4 << ','
          << r.Bx   << ',' << r.By << ',' << r.Bz << ','
          << r.mx   << ',' << r.my << ',' << r.mz
          << '\n';
    }
    csv.close();

    std::cout << "Wrote sim results to "
              << (fs::path(OUTPUT_DIR) / "sim_results.csv")
              << "\n";

    return 0;
}
