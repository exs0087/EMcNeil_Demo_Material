// include/Vec4.h
#pragma once

#include <cmath>
#include <cstddef>

/// 4‐component quaternion or general vector (scalar‐last: x,y,z,w)
struct Vec4 {
    double x{}, y{}, z{}, w{1.0};

    // Constructors
    constexpr Vec4() = default;
    constexpr Vec4(double _x, double _y, double _z, double _w)
      : x(_x), y(_y), z(_z), w(_w) {}

    // Element access
    double&       operator[](std::size_t i)       {
        switch(i){ case 0: return x; case 1: return y; case 2: return z; default: return w; }
    }
    double        operator[](std::size_t i) const {
        switch(i){ case 0: return x; case 1: return y; case 2: return z; default: return w; }
    }

    // Arithmetic
    Vec4 operator+(const Vec4& o) const { return {x+o.x, y+o.y, z+o.z, w+o.w}; }
    Vec4 operator*(double s)       const { return {x*s,   y*s,   z*s,   w*s  }; }
    friend Vec4 operator*(double s, const Vec4& v) { return v*s; }

    // Quaternion multiplication (Hamilton product)
    Vec4 operator*(const Vec4& o) const {
        return {
            w*o.x + x*o.w + y*o.z - z*o.y,
            w*o.y - x*o.z + y*o.w + z*o.x,
            w*o.z + x*o.y - y*o.x + z*o.w,
            w*o.w - x*o.x - y*o.y - z*o.z
        };
    }

    // Norm & normalize
    double norm() const {
        return std::sqrt(x*x + y*y + z*z + w*w);
    }
    void normalize() {
        double n = norm();
        if(n>0) { x/=n; y/=n; z/=n; w/=n; }
    }
    Vec4 normalized() const {
        Vec4 r = *this; r.normalize(); return r;
    }
};
