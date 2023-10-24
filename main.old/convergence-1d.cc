#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include "../quadrature/monte-carlo.h"
#include "../plot/convergence.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "output.svg";
	double max_samples_or_time = 1000;
	bool plot_samples = true;
	std::size_t seed = std::random_device()();
	auto [func,gt] = function1d(argc,argv);
	double rmin = 0; double rmax = 1;
	double min_plottable_error = -10.0;
    std::size_t montecarlo_error_average = 100;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-max-samples") { 
			plot_samples = true;
			max_samples_or_time = atof(argv[++i]); }
		else if (std::string(argv[i])=="-max-time") { 
			plot_samples = false;
			max_samples_or_time = atof(argv[++i]); } 			
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); }
	    else if ( (std::string(argv[i]) == "-range") && (i<(argc-2)) ) {
			rmin = atof(argv[++i]); rmax = atof(argv[++i]);
		}
		else if (std::string(argv[i])=="-min-plottable-error") { min_plottable_error = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-montecarlo-error-average") { montecarlo_error_average = atof(argv[++i]); } 
	}

	double ground_truth = gt(rmin,rmax);
	auto error = [ground_truth,min_plottable_error] (double a) {
		return std::max(min_plottable_error,std::abs(ground_truth-a));
	};

	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-60,-60,860,860));
	auto& graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({800,800},svg_cpp_plot::BoundingBox(0,std::log(min_plottable_error),double(max_samples_or_time),0)));
	graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.yticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xlabels(2);
	graph.ylabels(2);

    graph.area().add(convergence(plot_samples,stepper_monte_carlo_uniform(seed), func, range(rmin,rmax), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::red).stroke_width(2);
    graph.area().add(convergence(plot_samples,stepper_adaptive(nested(simpson,trapezoidal)), func, range(rmin, rmax), error, max_samples_or_time))
		.stroke(svg_cpp_plot::green).stroke_width(2);
    graph.area().add(convergence(plot_samples,stepper_adaptive(nested(boole,simpson)), func, range(rmin, rmax), error, max_samples_or_time))
		.stroke(svg_cpp_plot::blue).stroke_width(2);
	
    //Calculate function ranges
	float func_min = func(rmin); float func_max = func(rmin);
	for (float x = rmin; x <= rmax; x += ((rmax-rmin)/100.0f)) {
		if (func(x) < func_min) func_min = func(x);
		if (func(x) > func_max) func_max = func(x);
	}

	func_max += (func_max - func_min)*0.1f;
	func_min -= (func_max - func_min)*0.1f;

	auto& function_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(490,10))).add(svg_cpp_plot::Graph2D({300,200},svg_cpp_plot::BoundingBox(rmin,func_min,rmax,func_max)));
	function_graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	function_graph.area().add(svg_cpp_plot::_2d::function(func,rmin,rmax)).stroke(svg_cpp_plot::orange).stroke_width(1);

	
	std::ofstream f(output);
	f << svg;
	
	return 0;
}



