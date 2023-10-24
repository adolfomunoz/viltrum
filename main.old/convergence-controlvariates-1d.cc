#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include "../quadrature/monte-carlo.h"
#include "../quadrature/control-variates.h"
#include "../plot/convergence.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>


template<typename CV>
class FunctionCV {
    CV cv;
public:
    double operator()(double x) const {
        return cv(std::array<double,1>{x});
    }

    FunctionCV(CV&& c) : cv(std::forward<CV>(c)) { }
}; 

template<typename CV>
auto function_cv(CV&& cv) { return FunctionCV<CV>(std::forward<CV>(cv)); }

int main(int argc, char **argv) {	
	const char* output = "output.svg";
	unsigned long max_samples = 1000;
	std::size_t seed = std::random_device()();
	auto [func,gt] = function1d(argc,argv);
	double rmin = 0; double rmax = 1;
	double min_plottable_error = -10.0;
    std::size_t montecarlo_error_average = 100;
    unsigned long adaptive_iterations = 11;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-max-samples") { max_samples = atol(argv[++i]); } 
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); }
	        else if ( (std::string(argv[i]) == "-range") && (i<(argc-2)) ) {
			rmin = atof(argv[++i]); rmax = atof(argv[++i]);
		}
		else if (std::string(argv[i])=="-min-plottable-error") { min_plottable_error = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-montecarlo-error-average") { montecarlo_error_average = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-adaptive-iterations") { adaptive_iterations = atol(argv[++i]); } 
	}

	double ground_truth = gt(rmin,rmax);
	auto error = [ground_truth,min_plottable_error] (double a) {
		return std::max(min_plottable_error,std::abs(ground_truth-a));
	};

	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-60,-60,860,860));
	auto& graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({800,800},svg_cpp_plot::BoundingBox(0,std::log(min_plottable_error),double(max_samples),0)));
	graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.yticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xlabels(2);
	graph.ylabels(2);

    graph.area().add(convergence_samples(stepper_monte_carlo_uniform(seed), func, range(rmin,rmax), error, max_samples, montecarlo_error_average))
		.stroke(svg_cpp_plot::red).stroke_width(2);
    graph.area().add(convergence_samples(stepper_control_variate(control_variate_quadrature(trapezoidal),seed), func, range(rmin, rmax), error, max_samples, montecarlo_error_average))
		.stroke(svg_cpp_plot::green).stroke_width(2);
    graph.area().add(convergence_samples(stepper_control_variate(control_variate_quadrature(simpson),seed), func, range(rmin, rmax), error, max_samples, montecarlo_error_average))
		.stroke(svg_cpp_plot::blue).stroke_width(2);
    graph.area().add(convergence_samples(stepper_control_variate(control_variate_quadrature(boole),seed), func, range(rmin, rmax), error, max_samples, montecarlo_error_average))
		.stroke(svg_cpp_plot::purple).stroke_width(2);
    graph.area().add(convergence_samples(stepper_control_variate(control_variate_quadrature_adaptive(nested(simpson, trapezoidal),adaptive_iterations),seed), func, range(rmin, rmax), error, max_samples, montecarlo_error_average))
		.stroke(svg_cpp_plot::orange).stroke_width(2);

	
    //Calculate function ranges
	float func_min = func(rmin); float func_max = func(rmin);
	for (float x = rmin; x <= rmax; x += ((rmax-rmin)/100.0f)) {
		if (func(x) < func_min) func_min = func(x);
		if (func(x) > func_max) func_max = func(x);
	}

	func_max += (func_max - func_min)*0.1f;
	func_min -= (func_max - func_min)*0.1f;

    auto f = [&] (const std::array<double,1>& x) { return func(std::get<0>(x)); };
	auto& function_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(490,10))).add(svg_cpp_plot::Graph2D({300,200},svg_cpp_plot::BoundingBox(rmin,func_min,rmax,func_max)));
	function_graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	function_graph.area().add(svg_cpp_plot::_2d::function(function_cv(control_variate_quadrature(trapezoidal)(f,range(rmin,rmax))),rmin,rmax)).stroke(svg_cpp_plot::green).stroke_width(1);
	function_graph.area().add(svg_cpp_plot::_2d::function(function_cv(control_variate_quadrature(simpson)(f,range(rmin,rmax))),rmin,rmax)).stroke(svg_cpp_plot::blue).stroke_width(1);
	function_graph.area().add(svg_cpp_plot::_2d::function(function_cv(control_variate_quadrature(boole)(f,range(rmin,rmax))),rmin,rmax)).stroke(svg_cpp_plot::purple).stroke_width(1);
	function_graph.area().add(svg_cpp_plot::_2d::function(function_cv(control_variate_quadrature_adaptive(nested(simpson, trapezoidal),adaptive_iterations)(f,range(rmin,rmax))),rmin,rmax)).stroke(svg_cpp_plot::orange).stroke_width(1);
	function_graph.area().add(svg_cpp_plot::_2d::function(func,rmin,rmax)).stroke(svg_cpp_plot::red).stroke_width(1);

	
	std::ofstream saving_file(output);
	saving_file << svg;
	
	return 0;
}



