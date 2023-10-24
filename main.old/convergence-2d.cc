#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate.h"
#include "../quadrature/monte-carlo.h"
#include "../quadrature/control-variates.h"
#include "../plot/convergence.h"
#include "../plot/control-variate.h"
#include "../functions/functions2d.h"
#include "../quadrature/integrate-adaptive-control-variates.h"

#include <iostream>
#include <list>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "output.svg";
	double max_samples_or_time = 1000;
	bool plot_samples = true;
	std::size_t seed = std::random_device()();
	auto [func,gt] = function2d(argc,argv);
//	double rmin = 0; double rmax = 1;
	double min_plottable_error = -10.0;
	int plot_resolution = 300;
    std::size_t montecarlo_error_average = 100;
    unsigned long adaptive_iterations = 11;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-max-samples") { 
			plot_samples = true;
			max_samples_or_time = atof(argv[++i]); }
		else if (std::string(argv[i])=="-max-time") { 
			plot_samples = false;
			max_samples_or_time = atof(argv[++i]); } 			
		else if (std::string(argv[i])=="-seed") { seed = atol(argv[++i]); }
//	    else if ( (std::string(argv[i]) == "-range") && (i<(argc-2)) ) {
//			rmin = atof(argv[++i]); rmax = atof(argv[++i]);
//		}
		else if (std::string(argv[i])=="-min-plottable-error") { min_plottable_error = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-montecarlo-error-average") { montecarlo_error_average = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-plot-resolution") { plot_resolution = atoi(argv[++i]); }
		else if (std::string(argv[i])=="-adaptive-iterations") { adaptive_iterations = atol(argv[++i]); } 
	}

	double ground_truth = gt(0,0,1,1);
	auto error = [ground_truth,min_plottable_error] (double a) {
		return std::max(min_plottable_error,std::abs(ground_truth-a));
	};

	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-40,-40,1140,840));
	auto& graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({800,800},svg_cpp_plot::BoundingBox(0,std::log(min_plottable_error),max_samples_or_time,0)));
	graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.yticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xlabels(2);
	graph.ylabels(2);

    graph.area().add(convergence(plot_samples,stepper_monte_carlo_uniform(seed), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::red).stroke_width(2);
    graph.area().add(convergence(plot_samples,stepper_adaptive(nested(simpson,trapezoidal)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time))
		.stroke(svg_cpp_plot::rgb(0,0.5,0)).stroke_width(2);
    graph.area().add(convergence(plot_samples,stepper_adaptive(nested(boole,simpson)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time))
		.stroke(svg_cpp_plot::rgb(0,0,0.5)).stroke_width(2);
    graph.area().add(convergence(plot_samples,stepper_control_variate(control_variate_quadrature_adaptive(nested(simpson, trapezoidal),adaptive_iterations)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0,1,0)).stroke_width(2);
    graph.area().add(convergence(plot_samples,stepper_control_variate(control_variate_quadrature_adaptive(nested(boole, simpson),adaptive_iterations)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0,0,1)).stroke_width(2);
	graph.area().add(convergence(plot_samples,stepper_adaptive_control_variates(nested(simpson, trapezoidal), adaptive_iterations, seed, seed+1), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0,1,0)).stroke_width(2).stroke_dasharray({4,2});
	graph.area().add(convergence(plot_samples,stepper_adaptive_control_variates(nested(boole,simpson), adaptive_iterations, seed, seed+1), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0,0,1)).stroke_width(2).stroke_dasharray({4,2});
	
	auto& function_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(540,10))).add(svg_cpp_plot::Graph2D({250,250},svg_cpp_plot::BoundingBox(0,0,1,1)));
	function_graph.border().stroke(svg_cpp_plot::black).stroke_width(3);
    function_graph.area().add(svg_cpp_plot::_2d::function_2d(func,{0.0,0.0},{1.0,1.0},{plot_resolution,plot_resolution}));
	auto& function_mc_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(850,0))).add(svg_cpp_plot::Graph2D({250,250},svg_cpp_plot::BoundingBox(0,0,1,1)));
    function_mc_graph.border().stroke(svg_cpp_plot::red).stroke_width(3);
    function_mc_graph.area().add(svg_cpp_plot::_2d::function_2d([] (double x, double y) { return 0.0; },{0.0,0.0},{1.0,1.0},{plot_resolution,plot_resolution}));
	auto& function_cv_low_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(850,275))).add(svg_cpp_plot::Graph2D({250,250},svg_cpp_plot::BoundingBox(0,0,1,1)));
    function_cv_low_graph.border().stroke(svg_cpp_plot::rgb(0,1,0)).stroke_width(3);
    function_cv_low_graph.area().add(plot_control_variate(control_variate_quadrature_adaptive(nested(simpson, trapezoidal),adaptive_iterations),func,range(0.0,0.0,1.0,1.0),plot_resolution));
	auto& function_cv_high_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(850,550))).add(svg_cpp_plot::Graph2D({250,250},svg_cpp_plot::BoundingBox(0,0,1,1)));
    function_cv_high_graph.border().stroke(svg_cpp_plot::rgb(0,0,1)).stroke_width(3);
    function_cv_high_graph.area().add(plot_control_variate(control_variate_quadrature_adaptive(nested(boole, simpson),adaptive_iterations),func,range(0.0,0.0,1.0,1.0),plot_resolution));
                

	
	std::ofstream f(output);
	f << svg;
	
	return 0;
}



