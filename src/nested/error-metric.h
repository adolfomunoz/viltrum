#pragma once
#include <cmath>

namespace viltrum {

/**
 * Error metrics compare two elements and give an error. More of them can be developed, and you may need to develop your own ones for your own data types.
 */

class error_metric_absolute {
public:
    float operator()(float a, float b) const { return std::abs(b-a); }
    double operator()(double a, double b) const { return std::abs(b-a); }
    template<typename V>
    auto operator()(const V& v1, const V& v2, typename std::enable_if_t<std::is_arithmetic<V::value_type>, int> = 0) const {
        auto i1 = v1.begin(); auto i2 = v2.begin();
        using S = std::decay_t<decltype((*this)(*i1,*i2))>;
        S s; bool first = true;
        if ( (i1 != v1.end()) && (i2 != v2.end()) ) s = (*this)(*i1,*i2);
        while ( (i1 != v1.end()) && (i2 != v2.end()) ) {
            if (first) { s = (*this)(*i1,*i2); first = false; }
            else { s += (*this)(*i1,*i2); }
            ++i1; ++i2;
        }
        return s;
    }     
};

class error_metric_relative {
public:
    float operator()(float a, float b) const { return std::abs(b-a)/(std::max(std::abs(a),std::abs(b))); }
    double operator()(double a, double b) const { return std::abs(b-a)/(std::max(std::abs(a),std::abs(b))); }
    template<typename V>
    auto operator()(const V& v1, const V& v2, typename std::enable_if_t<std::is_arithmetic<V::value_type>, int> = 0) const {
        auto i1 = v1.begin(); auto i2 = v2.begin();
        using S = std::decay_t<decltype((*this)(*i1,*i2))>;
        S s; bool first = true;
        if ( (i1 != v1.end()) && (i2 != v2.end()) ) s = (*this)(*i1,*i2);
        while ( (i1 != v1.end()) && (i2 != v2.end()) ) {
            if (first) { s = (*this)(*i1,*i2); first = false; }
            else { s += (*this)(*i1,*i2); }
            ++i1; ++i2;
        }
        return s;
    }     
};

};