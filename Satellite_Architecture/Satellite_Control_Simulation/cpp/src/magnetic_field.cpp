#include "magnetic_field.hpp"

// TODO: replace these stubs with your MATLAB‐ported code

Vec3 issOrbitECI(double /*t*/) {
    // stub: circular at 7000 km altitude along x‐axis
    return {7000e3, 0.0, 0.0};
}

Vec3 igrfECEF(double /*t*/, double /*lat*/, double /*lon*/, double /*alt*/) {
    // stub: constant field vector [nT]
    return {20000.0e-9, 0.0, 0.0};  // [T]
}
