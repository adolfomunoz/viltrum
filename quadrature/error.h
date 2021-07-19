#pragma once

#include <cmath>

namespace viltrum {

class ErrorSingleDimensionStandard {
    float norm(float f) const { return std::abs(f); }
    double norm(double f) const { return std::abs(f); }
    template<typename V>
    auto norm(const V& v) const {
        auto i = v.begin();
        auto s = norm(*i); ++i;
        while (i!=v.end()) {
            s += norm(*i); ++i;
        }
        return s;
    }

public:
    template<typename R>
    auto operator()(const R& region) const {
        auto max_err = norm(region.error(0)); std::size_t max_dim = 0; 
        decltype(max_err) err;
        for (std::size_t d = 1; d<R::dimensions; ++d) {
            err = norm(region.error(d));
            if (err>max_err) {
                max_err = err; max_dim = d;
            }
        }
        return std::make_tuple(max_err,max_dim);
    }
};

ErrorSingleDimensionStandard error_single_dimension_standard() { return ErrorSingleDimensionStandard(); }
ErrorSingleDimensionStandard error_absolute_single_dimension() { return ErrorSingleDimensionStandard(); }

class ErrorRelativeSingleDimensionStandard {
    float norm(float f) const { return std::abs(f); }
    double norm(double f) const { return std::abs(f); }
    template<typename V>
    auto norm(const V& v) const {
        auto i = v.begin();
        auto s = norm(*i); ++i;
        while (i!=v.end()) {
            s += norm(*i); ++i;
        }
        return s;
    }
    
    double offset;

public:
    template<typename R>
    auto operator()(const R& region) const {
        using real = decltype(norm(region.integral()));
        double den = std::max(norm(region.integral()),real(offset));
        auto max_err = norm(region.error(0))/den; std::size_t max_dim = 0; 

        decltype(max_err) err;
        for (std::size_t d = 1; d<R::dimensions; ++d) {
            err = norm(region.error(d))/den;
            if (err>max_err) {
                max_err = err; max_dim = d;
            }
        }
        return std::make_tuple(max_err,max_dim);
    }
    
    ErrorRelativeSingleDimensionStandard(double offset) : offset(offset) {}
};

ErrorRelativeSingleDimensionStandard error_relative_single_dimension(double offset = 1.e-6) { return ErrorRelativeSingleDimensionStandard(offset); }


class ErrorSingleDimensionSize {
    float norm(float f) const { return std::abs(f); }
    double norm(double f) const { return std::abs(f); }
    template<typename V>
    auto norm(const V& v) const {
        auto i = v.begin();
        auto s = norm(*i); ++i;
        while (i!=v.end()) {
            s += norm(*i); ++i;
        }
        return s;
    }

    double size_factor;

public:
    template<typename R>
    auto operator()(const R& region) const {
        auto max_err = norm(region.error(0)); std::size_t max_dim = 0;
        max_err += std::abs(size_factor*(region.range().max(0) - region.range().min(0))); 
        decltype(max_err) err;
        for (std::size_t d = 1; d<R::dimensions; ++d) {
            err = norm(region.error(d));
            err += std::abs(size_factor*(region.range().max(d) - region.range().min(d)));
            if (err>max_err) {
                max_err = err; max_dim = d;
            }
        }
        return std::make_tuple(max_err,max_dim);
    }

    ErrorSingleDimensionSize(double sf) : size_factor(sf) { }
};

ErrorSingleDimensionSize error_single_dimension_size(double size_factor = 0.01) { return ErrorSingleDimensionSize(size_factor); }

class ErrorRelativeSingleDimensionSize {
    float norm(float f) const { return std::abs(f); }
    double norm(double f) const { return std::abs(f); }
    template<typename V>
    auto norm(const V& v) const {
        auto i = v.begin();
        auto s = norm(*i); ++i;
        while (i!=v.end()) {
            s += norm(*i); ++i;
        }
        return s;
    }

    double size_factor;
    double offset;

public:
    template<typename R>
    auto operator()(const R& region) const {
        double den = std::max(norm(region.integral()),offset);
        auto max_err = norm(region.error(0))/den; std::size_t max_dim = 0;
        max_err += std::abs(size_factor*(region.range().max(0) - region.range().min(0)));
        decltype(max_err) err;
        for (std::size_t d = 1; d<R::dimensions; ++d) {
            err = norm(region.error(d))/den;
            err += std::abs(size_factor*(region.range().max(d) - region.range().min(d)));
            if (err>max_err) {
                max_err = err; max_dim = d;
            }
        }
        return std::make_tuple(max_err,max_dim);
    }

    ErrorRelativeSingleDimensionSize(double sf, double offset) : size_factor(sf),offset(offset) { }
};

ErrorRelativeSingleDimensionSize error_relative_single_dimension_size(double size_factor = 1.e-5, double offset = 1.e-6) { return ErrorRelativeSingleDimensionSize(size_factor, offset); }

class ErrorPartiallyRelativeSingleDimensionSize {
    float norm(float f) const { return std::abs(f); }
    double norm(double f) const { return std::abs(f); }
    template<typename V>
    auto norm(const V& v) const {
        auto i = v.begin();
        auto s = norm(*i); ++i;
        while (i!=v.end()) {
            s += norm(*i); ++i;
        }
        return s;
    }

    double size_factor;
    std::size_t relative_dimensions;
    double offset;

public:
    template<typename R>
    auto operator()(const R& region) const {
        double den = std::max(norm(region.integral()),offset);
        auto max_err = norm(region.error(0)); std::size_t max_dim = 0;
        if (relative_dimensions>0) max_err/=den;
        max_err += std::abs(size_factor*(region.range().max(0) - region.range().min(0)));
        decltype(max_err) err;
        for (std::size_t d = 1; d<R::dimensions; ++d) {
            err = norm(region.error(d));
            if (d<relative_dimensions) err/=den;
            err += std::abs(size_factor*(region.range().max(d) - region.range().min(d)));
            if (err>max_err) {
                max_err = err; max_dim = d;
            }
        }
        return std::make_tuple(max_err,max_dim);
    }

    ErrorPartiallyRelativeSingleDimensionSize(double sf, std::size_t relative_dimensions, double offset) : size_factor(sf),relative_dimensions(relative_dimensions),offset(offset) { }
};

ErrorPartiallyRelativeSingleDimensionSize error_partially_relative_single_dimension_size(double size_factor = 1.e-5, std::size_t relative_dimensions = 2, double offset = 1.e-6) { return ErrorPartiallyRelativeSingleDimensionSize(size_factor, relative_dimensions, offset); }

}





