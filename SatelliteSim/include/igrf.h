// include/igrf.h
#pragma once

#include <string>
#include <vector>
#include "Vec3.h"

void geodeticToGeocentric(double lat_rad,
                          double lon_rad,
                          double alt_km,
                          double &r_km,
                          double &cos_theta,
                          double &sin_theta);

std::tuple<double,double,double>
computeSphericalField(double r_km,
                      double theta_rad,
                      double phi_rad,
                      const std::vector<double> &gh);

Vec3 igrf(double time,
          double lat_rad,
          double lon_rad,
          double alt_km,
          const std::string &coord = "geodetic");
