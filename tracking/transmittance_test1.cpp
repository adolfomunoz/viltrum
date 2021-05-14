#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../tracking/ratio-tracking.h"
#include "../tracking/delta_tracking.h"
#include "../tracking/ray_marching.h"
#include "../quadrature/integrate.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

typedef Range<float,1> RangeF;

int main(int argc, char **argv) 
{	
	bool verbose = false;
	unsigned int max_log2_samples = 10;
	int seed = 1;

	for (int i = 0; i < argc; ++i) {
		if (std::string(argv[i]) == "-verbose")
			verbose = true;	

		if (std::string(argv[i]) == "-max_log2_samples" && i < argc-1)
			max_log2_samples = atoi(argv[i+1]);		

		if (std::string(argv[i]) == "-seed" && i < argc-1)
			seed = atoi(argv[i+1]);		
	}


	auto [func,gt] = function1d(argc,argv);


	for(unsigned int i = 1; i<max_log2_samples; ++i)
	{

		unsigned long sp = pow(2,i);


		// Ray Marching
		auto v_raymarching(ray_marching(sp));
		using ResultRM = decltype(v_raymarching.init(func, range(0.f,1.f)));
		ResultRM v_rmsamples;

		//const F& f, const Range<Float,1>& range, Samples<Result>& samples
		v_raymarching.step(func, range(0.f,1.f), v_rmsamples, verbose);
		float rmintegral = v_raymarching.integral(func, v_rmsamples);
		unsigned long rmnb_queries = v_raymarching.queries(func, v_rmsamples);
		printf("RayMarching\t- Integral: %.10f [%d queries] \n", rmintegral, rmnb_queries);


		// Ratio Tracking
		auto v_ratiotracking(ratio_tracking(std::mt19937_64(seed),2.,0.));

		using ResultRatTracking = decltype(v_ratiotracking.init(func, range(0.f,1.f)));
		ResultRatTracking v_rtsamples;

		for(unsigned long j=0; j<sp; ++j)
			v_ratiotracking.step(func, range(0.f,1.f), v_rtsamples, verbose);
	
		float rtintegral = v_ratiotracking.integral(func, v_rtsamples);
		unsigned long rtnb_queries = v_ratiotracking.queries(func, v_rtsamples);


		printf("RatioTracking\t- Integral: %.10f [%d samples - %d queries] \n", rtintegral, sp, rtnb_queries);


		// Delta Tracking
		auto v_deltatracking(delta_tracking(std::mt19937_64(seed),2.));

		using ResultDeltaTracking = decltype(v_deltatracking.init(func, range(0.f,1.f)));
		ResultDeltaTracking v_dtsamples;

		for(unsigned long j=0; j<sp; ++j)
			v_deltatracking.step(func, range(0.f,1.f), v_dtsamples, verbose);
	
		float dtintegral = v_deltatracking.integral(func, v_dtsamples);
		unsigned long dtnb_queries = v_deltatracking.queries(func, v_dtsamples);


		printf("DeltaTracking\t- Integral: %.10f [%d samples - %d queries] \n", dtintegral, sp, dtnb_queries);

	}
	/*const char* output = "output.svg";
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
		return std::max(min_plottable_error,std::log(std::abs(ground_truth-a)));
	};

	
	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-60,-60,860,860));
	auto& graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({800,800},svg_cpp_plot::BoundingBox(0,min_plottable_error,double(max_samples),0)));
	graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.yticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xlabels(2);
	graph.ylabels(2);

    graph.area().add(convergence_samples(RayMarching(.01), func, range(rmin,rmax), error, max_samples, montecarlo_error_average))
		.stroke(svg_cpp_plot::red).stroke_width(2);
	*/


/*    graph.area().add(convergence_samples(stepper_control_variate(control_variate_quadrature(trapezoidal),seed), func, range(rmin, rmax), error, max_samples, montecarlo_error_average))
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
	function_graph.area().add(svg_cpp_plot::_2d::function(func,rmin,rmax)).stroke(svg_cpp_plot::red).stroke_width(1);*/

	/*
	std::ofstream saving_file(output);
	saving_file << svg;
	
	return 0;*/
}



