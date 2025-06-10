import os, sys

# Add the C++ build folder to Python’s path
root = os.path.dirname(__file__)
build_dir = os.path.abspath(os.path.join(root, '..', 'cpp', 'cmake-build-debug-wsl'))
sys.path.insert(0, build_dir)

import numpy as np
import satellite_kernels as sk

# Example: ECEF -> ENU -> spherical
r_ecef = np.array([3.0, 2.0, 1.0])
lat, lon = 38.0, 60.0

r_enu = sk.enu_from_ecef(r_ecef, lat, lon)
print("ENU vector:", r_enu)

range_, az, el = sk.enu_to_spherical(r_enu)
print(f"Range={range_:.2f} m, Az={az:.2f}°, El={el:.2f}°")

# Example: LST
jd_example = 2459930.5  # Example Julian Date
lon_obs = -104.0        # Observer longitude (deg)
lst = sk.lst_from_jd_lon(jd_example, lon_obs)
print(f"Local Sidereal Time: {lst:.2f}°")
