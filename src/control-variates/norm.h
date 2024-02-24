#pragma once
#include <cmath>
#include <type_traits>



namespace viltrum {
struct NormDefault {
    float operator()(float v) const { return std::abs(v); }
    double operator()(double v) const { return std::abs(v); }
    template<typename V>
    auto operator()(const V& v, typename std::enable_if_t<std::is_arithmetic_v<typename V::value_type>, int> = 0) const {
        auto i = v.begin(); 
        using S = decltype(norm(*i));
        S s; bool first = true;
        if (i != v.end()) { s = norm(*i); ++i; }
        while (i != v.end()) { s += norm(*i); ++i; }
        return s;
    }     
};
} 
