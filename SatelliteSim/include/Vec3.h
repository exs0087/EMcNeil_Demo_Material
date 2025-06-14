// include/Vec3.h
#pragma once

#include <cmath>
#include <cstddef>

/// 3‐component Cartesian vector
struct Vec3 {
    double x{}, y{}, z{};

    // Constructors
    constexpr Vec3() = default;
    constexpr Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

    // Element access
    double&       operator[](std::size_t i)       { return i==0?x : i==1?y : z; }
    double        operator[](std::size_t i) const { return i==0?x : i==1?y : z; }

    // Arithmetic
    Vec3 operator+(const Vec3& o)  const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o)  const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s)        const { return {x*s,   y*s,   z*s  }; }
    Vec3 operator/(double s)        const { return {x/s,   y/s,   z/s  }; }

    Vec3& operator+=(const Vec3& o)      { x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vec3& operator-=(const Vec3& o)      { x-=o.x; y-=o.y; z-=o.z; return *this; }
    Vec3& operator*=(double s)           { x*=s;   y*=s;   z*=s;   return *this; }
    Vec3& operator/=(double s)           { x/=s;   y/=s;   z/=s;   return *this; }

    // Dot & cross
    double dot(const Vec3& o)    const { return x*o.x + y*o.y + z*o.z; }
    Vec3   cross(const Vec3& o)  const {
        return { y*o.z - z*o.y,
                z*o.x - x*o.z,
                x*o.y - y*o.x };
    }

    // Norms
    double norm()        const { return std::sqrt(dot(*this)); }
    Vec3   normalized()  const {
        double n = norm();
        return n>0 ? (*this)/n : Vec3{};
    }
};
