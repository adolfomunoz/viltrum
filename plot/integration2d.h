#pragma once

#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include <cmath>


template<typename F,typename P>
class FunctionWrapper2d {
	F f;
	P& points;
	mutable std::shared_ptr<unsigned long> nsamples;
public:
	FunctionWrapper2d(const F& f, P& points) : f(f), points(points), 
		nsamples(std::make_shared<unsigned long>(0)) {}
	FunctionWrapper2d(F&& f, P& points) : f(std::forward<F>(f)), points(points), 
		nsamples(std::make_shared<unsigned long>(0)) {}
	
	auto operator()(const std::array<double,2>& pos) const {
		auto value = f(pos);
		points.add_point(std::get<0>(pos),std::get<1>(pos));
		++(*nsamples);
		return value;		
	}
	
	auto operator()(const std::array<float,2>& pos) const {
		auto value = f(pos);
		points.add_point(std::get<0>(pos),std::get<1>(pos));
		++(*nsamples);
		return value;		
	}
		
	unsigned long samples() const { return *nsamples; }
};

template<typename IntegratorBins, typename F, typename Float>
auto plot_samples_2d(const IntegratorBins& integrator, const F& function, const Range<Float,2>& range, unsigned long bins = 1000) {
	auto points = svg_cpp_plot::_2d::points();
	FunctionWrapper2d f(function,points);
	std::vector<decltype(f(range.min()))> result(bins,0);
    auto vb = vector_bins(result);
    integrator.integrate(vb, vector_resolution(result), f, range);
	return points;
}

template<typename Nested, typename Error, typename F, typename Float>
svg_cpp_plot::_2d::Group plot_adaptive_boundaries_2d(Nested&& nested, Error&& error, const F& f, const Range<Float,2>& range, unsigned long iterations) {
	auto cv_stepper = stepper_adaptive(std::forward<Nested>(nested),std::forward<Error>(error));
	auto regions = cv_stepper.init(f,range);
    for (unsigned long i = 0; i<iterations; ++i) cv_stepper.step(f,range,regions);
	auto sol = svg_cpp_plot::_2d::group();
	for (const auto& r : regions)
		sol.add(svg_cpp_plot::_2d::rect(
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)}));
	return sol;
}

template<typename Nested, typename Error, typename F, typename Float>
auto plot_adaptive_samples_2d(Nested&& nested, Error&& error, const F& function, const Range<Float,2>& range, unsigned long iterations) {
	auto points = svg_cpp_plot::_2d::points();
	FunctionWrapper2d f(function,points);
	auto cv_stepper = stepper_adaptive(std::forward<Nested>(nested),std::forward<Error>(error));
	auto regions = cv_stepper.init(f,range);
    for (unsigned long i = 0; i<iterations; ++i) cv_stepper.step(f,range,regions);
	return points;
}

template<typename Nested, typename Error, typename F, typename Float, typename ColorMap>
svg_cpp_plot::_2d::Group plot_adaptive_approximation_2d(Nested&& nested, Error&& error, const F& f, const Range<Float,2>& range, unsigned long iterations, const ColorMap& color_map, int resolution = 200) {
	auto cv_stepper = stepper_adaptive(std::forward<Nested>(nested),std::forward<Error>(error));
	auto regions = cv_stepper.init(f,range);
    for (unsigned long i = 0; i<iterations; ++i) cv_stepper.step(f,range,regions);
	auto sol = svg_cpp_plot::_2d::group();
	for (const auto& r : regions) 
		sol.add(svg_cpp_plot::_2d::function_2d(
			[&r] (Float x, Float y) {  return r.approximation_at(std::array{x,y}); },
			color_map,
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)},
			{int(resolution*(r.range().max(0)-r.range().min(0))),
			int(resolution*(r.range().max(1)- r.range().min(1)))}));
	return sol;
}

namespace {

	inline std::tuple<float,float,float> to_array(const std::tuple<float,float,float>& t) { return t; }
	
	#ifdef EIGEN_CORE_H
	template<typename D>
	std::tuple<float,float,float> to_array(const Eigen::ArrayBase<D>& t) { 
		return std::tuple<float,float,float>{t.coeff(0),t.coeff(1),t.coeff(2)}; 
	}
	#endif
}

template<typename Nested, typename Error, typename F, typename Float>
svg_cpp_plot::_2d::Group plot_adaptive_approximation_2d(Nested&& nested, Error&& error, const F& f, const Range<Float,2>& range, unsigned long iterations, int resolution = 200) {
	auto cv_stepper = stepper_adaptive(std::forward<Nested>(nested),std::forward<Error>(error));
	auto regions = cv_stepper.init(f,range);
    for (unsigned long i = 0; i<iterations; ++i) cv_stepper.step(f,range,regions);
	auto sol = svg_cpp_plot::_2d::group();
	for (const auto& r : regions) 
		sol.add(svg_cpp_plot::_2d::function_image(
			[&r] (Float x, Float y) {  return to_array(r.approximation_at(std::array{x,y})); },
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)},
			{int(resolution*(r.range().max(0)-r.range().min(0))),
			int(resolution*(r.range().max(1)- r.range().min(1)))}));
	return sol;
}

template<typename Nested, typename Error, typename F, typename Float, typename ColorMap>
svg_cpp_plot::_2d::Group plot_adaptive_error_2d(Nested&& nested, Error&& error, const F& f, const Range<Float,2>& range, unsigned long iterations, const ColorMap& color_map, int resolution = 200) {
	auto cv_stepper = stepper_adaptive(std::forward<Nested>(nested),std::forward<Error>(error));
	auto regions = cv_stepper.init(f,range);
    for (unsigned long i = 0; i<iterations; ++i) cv_stepper.step(f,range,regions);
	auto sol = svg_cpp_plot::_2d::group();
	for (auto r : regions)
		sol.add(svg_cpp_plot::_2d::function_2d(
			[&f,&r] (Float x, Float y) {  return r.approximation_at(std::array{x,y}) - f(std::array{x,y}); },
			color_map,
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)},
			{int(resolution*(r.range().max(0)-r.range().min(0))),
			int(resolution*(r.range().max(1)- r.range().min(1)))}));
	return sol;
}




