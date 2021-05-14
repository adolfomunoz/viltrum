#pragma once

#include <svg-cpp-plot/svg-cpp-plot.h>
#include <vector>


template<typename CVG, typename F, typename Float>
auto plot_control_variate(const CVG& cv_generator, const F& f, const Range<Float,2>& range, int resolution = 200) {
    auto cv = cv_generator([f] (const std::array<Float,2>& x) { return f(x[0],x[1]); },range);
    return svg_cpp_plot::_2d::function_2d([cv] (Float x, Float y) { return cv(std::array<Float,2>{x,y}); },{range.min(0),range.min(1)},{range.max(0),range.max(1)},{resolution,resolution});
}


