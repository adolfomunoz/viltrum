#pragma once

#include <svg-cpp-plot/svg-cpp-plot.h>
#include <vector>

template<typename R>
svg_cpp_plot::_2d::Group plot_regions_2d(const std::vector<R>& regions, int resolution = 200) {
	auto sol = svg_cpp_plot::_2d::group();;
	auto& graph = sol.add(svg_cpp_plot::_2d::group());
	auto& boundaries = sol.add(svg_cpp_plot::_2d::group());
	for (auto r : regions) {
		graph.add(svg_cpp_plot::_2d::function_2d(
			[&r] (double x, double y) { return r.approximation_at(std::array{x,y}); },
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)},
			{int(resolution*abs(r.range().max(0)-r.range().min(0))),
			int(resolution*abs(r.range().max(1)- r.range().min(1)))}));
		boundaries.add(svg_cpp_plot::_2d::rect(
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)})).stroke(svg_cpp_plot::green).stroke_width(1);
	}
	return sol;
}

template<typename F, typename R>
svg_cpp_plot::_2d::Group plot_regions_error_2d(const F& f, const std::vector<R>& regions, int resolution = 200) {
	auto sol = svg_cpp_plot::_2d::group();
	for (auto r : regions) {
		sol.add(svg_cpp_plot::_2d::function_2d(
			[&r,&f] (double x, double y) { return r.approximation_at(std::array{x,y}) - f(x,y); }, svg_cpp_plot::_2d::color_map_red_blue(-1,1),
			{r.range().min(0),r.range().min(1)},
			{r.range().max(0),r.range().max(1)},
			{int(resolution*abs(r.range().max(0)-r.range().min(0))),
			int(resolution*abs(r.range().max(1)- r.range().min(1)))}));
	}
	return sol;
}
