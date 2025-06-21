// src/main.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <filesystem>
#include <string>
#include <cmath>

#include "adcs_hal.h"         // HAL_init, HAL_set_control_mode, ControlMode
#include "Vec3.h"
#include "Vec4.h"
#include "Vec7.h"
#include "eulerseqns2.h"      // eulerseqns2(), clearActHistory(), getActHistory()
#include "ode4.h"
#include "magnetic_field.h"   // getBodyFieldAt(time, Vec4{lat,lon,alt_km,0})
#include "iss_orbit.h"        // issOrbitECI(double t)

#ifndef OUTPUT_DIR
# define OUTPUT_DIR "./Test_Setup/cpp"
#endif
static constexpr char const* kOutputDir = OUTPUT_DIR;

// IGRF reference radius [km]
static constexpr double Re_km = 6371.2;

struct Row {
    double time;
    double wx, wy, wz;
    double q1, q2, q3, q4;
    double Bx, By, Bz;
    double mx, my, mz;
};

int main(int argc, char** argv) {
    // 1 Hz control update
    const double sampleRate = 1.0;    // seconds
    const double t_final   = 16000.0; // seconds
    const size_t N = static_cast<size_t>(t_final / sampleRate) + 1;

    // 1) Initialize HAL & clear any previous history
    HAL_init(sampleRate);
    clearActHistory();

    // 2) Choose control law
    if (argc > 1 && std::string(argv[1]) == "bdot") {
        HAL_set_control_mode(CTRL_MODE_B_DOT);
        std::cout << "Starting in Ḃ‐law (bang–bang) mode\n";
    } else {
        HAL_set_control_mode(CTRL_MODE_WX_B);
        std::cout << "Starting in ω×B proportional mode\n";
    }

    // 3) Time vector at 1 Hz
    std::vector<double> times(N);
    for (size_t i = 0; i < N; ++i) times[i] = i * sampleRate;

    // 4) Initial state quaternion + ang. rates
    Vec7 y0{0,0.2,0.2, 0,0,0,1};

    // 5) Integrate dynamics (RK4) — one control call per second
    auto t_start = std::chrono::high_resolution_clock::now();
    auto Y       = ode4(eulerseqns2, times, y0);
    auto t_end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << "Simulation took " << elapsed.count() << " s\n";

    // 6) Grab commanded dipole history
    auto const& act_hist = getActHistory();

    // 7) Assemble results
    std::vector<Row> allResults;
    allResults.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        Row r;
        r.time = times[i];

        // attitude states
        r.wx = Y[i][0]; r.wy = Y[i][1]; r.wz = Y[i][2];
        r.q1 = Y[i][3]; r.q2 = Y[i][4]; r.q3 = Y[i][5]; r.q4 = Y[i][6];

        // --- 8) get ECI position [km] for orbit propagator ---
        Vec3 posEci = issOrbitECI(times[i]);

        // 9) convert ECI → geodetic lat/lon/alt_km
        double x   = posEci.x;
        double y   = posEci.y;
        double z   = posEci.z;
        double rkm = std::hypot(x, y, z);
        double lon = std::atan2(y, x);           // radians
        double lat = std::asin(z / rkm);         // geocentric initial guess
        // ellipsoid conversion inside getBodyFieldAt now does full geodetic

        // 10) query IGRF with geodetic lat/lon/alt_km
        Vec3 Beci = getBodyFieldAt(
            times[i],
            Vec4{ lat, lon, rkm - Re_km, 0.0 }
        );
        r.Bx = Beci.x;
        r.By = Beci.y;
        r.Bz = Beci.z;

        // 11) dipole commands
        r.mx = act_hist[i].mx;
        r.my = act_hist[i].my;
        r.mz = act_hist[i].mz;

        allResults.push_back(r);
    }

    // 12) Write CSV
    namespace fs = std::filesystem;
    fs::create_directories(kOutputDir);
    fs::path out = fs::path(kOutputDir) / "sim_results.csv";
    std::ofstream csv{ out };
    csv << "time,wx,wy,wz,q1,q2,q3,q4,Bx,By,Bz,mx,my,mz\n";
    for (auto const& r : allResults) {
        csv
          << r.time << ','
          << r.wx   << ',' << r.wy << ',' << r.wz << ','
          << r.q1   << ',' << r.q2 << ',' << r.q3 << ',' << r.q4 << ','
          << r.Bx   << ',' << r.By << ',' << r.Bz << ','
          << r.mx   << ',' << r.my << ',' << r.mz << '\n';
    }
    csv.close();

    std::cout << "Wrote sim results to " << out << "\n";
    return 0;
}
