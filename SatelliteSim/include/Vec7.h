// include/Vec7.h
#pragma once

#include <cstddef>
#include <initializer_list>
#include <algorithm>

/// 7‐component general vector ([ω₁,ω₂,ω₃,q₁,q₂,q₃,q₄])
struct Vec7 {
    double v[7];

    // Zero‐initialize
    Vec7() { std::fill(v, v+7, 0.0); }

    // Init‐list constructor (pads with zeros if fewer than 7)
    Vec7(std::initializer_list<double> list) {
        std::size_t i=0;
        for(double x : list) {
            if(i<7) v[i++] = x;
        }
        for(; i<7; ++i) v[i] = 0.0;
    }

    // Element access
    double&       operator[](std::size_t i)       { return v[i]; }
    double        operator[](std::size_t i) const { return v[i]; }

    // Arithmetic
    Vec7 operator+(const Vec7& o) const {
        Vec7 r;
        for(std::size_t i=0;i<7;++i) r.v[i] = v[i] + o.v[i];
        return r;
    }
    Vec7 operator*(double s) const {
        Vec7 r;
        for(std::size_t i=0;i<7;++i) r.v[i] = v[i] * s;
        return r;
    }
    friend Vec7 operator*(double s, const Vec7& a) {
        return a * s;
    }
};
