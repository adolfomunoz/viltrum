#pragma once
#include <tuple>
#include <array>
#include <list>
#include <vector>

namespace viltrum {
 
/**
 * This function wrapper adapts functions of "N" parameters into the "array" version that is required for viltrum ndimensional integrators
 **/
template<typename F>
class FunctionWrapper {
	F f;
public:
	FunctionWrapper(const F& f) : f(f) {}
	FunctionWrapper(F&& f) : f(std::forward<F>(f)) {}
	
    template<typename Float, std::size_t N>
	auto operator()(const std::array<Float, N>& params) const {
		return std::apply(f,params);
	}
};

template<typename F>
auto function_wrapper(F&& f) {
    return FunctionWrapper<std::decay_t<F>>(std::forward<std::decay_t<F>>(f)); 
}

template<typename F>
auto function_wrapper(const F& f) {
    return FunctionWrapper<std::decay_t<F>>(f); 
}

/*
template<typename F>
auto function_wrapper(F& f) {
    return FunctionWrapper<std::decay_t<F>&>(f); 
}
*/





/**
 * This function wrapper helps counting number of function invocations. 
 * It also adapts functions of "N" parameters into the "array" version that is required for viltrum ndimensional integrators
 **/
template<typename F>
class FunctionWrapperCount {
	F f;
	mutable unsigned long evals = 0;
public:
	FunctionWrapperCount(const F& f) : f(f) {}
	FunctionWrapperCount(F&& f) : f(std::forward<F>(f)) {}
	
    template<typename Float, std::size_t N>
	auto operator()(const std::array<Float, N>& params) const {
		++evals;
		return std::apply(f,params);
	}
	
    unsigned long evaluations() const { return evals; }
};

/**
 * This function wrapper helps storing function invocations for plotting samples (for instance).
 * It also adapts functions of "N" parameters into the "array" version that is required for viltrum ndimensional integrators
 **/
template<typename F>
class FunctionWrapperProfile {
	F f;
    mutable std::vector<std::list<float>> params_;
    mutable std::list<float> values_;
public:
	FunctionWrapperProfile(const F& f) : f(f) {}
	FunctionWrapperProfile(F&& f) : f(std::forward<F>(f)) {}
	
    template<typename Float, std::size_t N>
	auto operator()(const std::array<Float, N>& params) const {
        if (N>params_.size()) params_.resize(N);
        for (std::size_t i = 0; i<N; ++i) params_[i].push_back(params[i]);
        auto sol = std::apply(f,params);
        values_.push_back(float(sol));
		return sol;
	}
	
    const std::list<float>& values() const { return values_; }
    const std::list<float>& params(std::size_t i) const { return params_[i%params_.size()]; }
    unsigned long evaluations() const { return values_.size(); }
};

};
