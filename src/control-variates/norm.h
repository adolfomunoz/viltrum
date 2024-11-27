#pragma once
#include <cmath>
#include <type_traits>
#include <complex>





namespace viltrum {
struct NormDefault {
    float operator()(float v) const { return std::abs(v); }
    double operator()(double v) const { return std::abs(v); }
    template<typename N>
    auto operator()(const std::complex<N>& c) const {
        return (*this)(c.real()) + (*this)(c.imag());
    } 
    template<typename V>
    auto operator()(const V& v) const {
        auto i = v.begin(); 
        using S = decltype((*this)(*i));
        S s; 
        if (i != v.end()) { s = (*this)(*i); ++i; }
        while (i != v.end()) { s += (*this)(*i); ++i; }
        return s;
    }
    //These sign functions are for when the sign matters (intermediate calculations) 
    float sign(float v) const { return v; }      
    float sign(double v) const { return v; }      
    template<typename N>
    auto sign(const std::complex<N>& c) const {
        return sign(c.real()) + sign(c.imag());
    } 
    template<typename V>
    auto sign(const V& v) const {
        auto i = v.begin(); 
        using S = decltype(sign(*i));
        S s; 
        if (i != v.end()) { s = sign(*i); ++i; }
        while (i != v.end()) { s += sign(*i); ++i; }
        return s;
    }};
} 
