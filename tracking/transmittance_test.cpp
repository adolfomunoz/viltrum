#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>

#include "../tracking/residual_ratio_tracking.h"
#include "../tracking/delta_tracking.h"
#include "../tracking/weighted_delta_tracking.h"
#include "../tracking/ray_marching.h"
#include "../quadrature/integrate.h"
#include "../functions/functions1d.h"

#include <iostream>
#include <cmath>

typedef Range<float,1> RangeF;

enum TrackingMethod {	TMRayMarching = 1, TMDeltaTracking = 2, TMWeightedDeltaTracking = 4, TMRatioTracking = 8, 
						TMResidualRatioTracking = 16, TMQuadrature = 32, TMCVQuadrature = 64,  TMAll = 255};

template<typename F,std::size_t DIM>
class FunctionWrapper {
	F f;
	mutable std::shared_ptr<unsigned long> nsamples;
public:
	FunctionWrapper(const F& f) : 
		f(f), nsamples(std::make_shared<unsigned long>(0)) {}
	FunctionWrapper(F&& f) : 
		f(std::forward<F>(f)), nsamples(std::make_shared<unsigned long>(0)) {}
	
	double operator()(const std::array<double,DIM>& pos) const {
		double value = std::apply(f,pos);
		++(*nsamples);
		return value;
	}
	
	unsigned long samples() const { return *nsamples; }
};



int main(int argc, char **argv) 
{	
	bool verbose = false;
	unsigned int max_log2_samples = 10;
	int seed = 1;

	const char* output = "output.svg";

	char method = 0;

	for (int i = 0; i < argc; ++i) {
		if (std::string(argv[i]) == "-verbose")
			verbose = true;	

		if (std::string(argv[i]) == "-max_log2_samples" && i < argc-1)
			max_log2_samples = atoi(argv[i+1]);		

		if (std::string(argv[i]) == "-seed" && i < argc-1)
			seed = atoi(argv[i+1]);	

		if (std::string(argv[i]) == "-raymarching")
			method = method | TrackingMethod::TMRayMarching; 	

		if (std::string(argv[i]) == "-deltatracking")
			method = method | TrackingMethod::TMDeltaTracking; 	

		if (std::string(argv[i]) == "-weighteddeltatracking")
			method = method | TrackingMethod::TMWeightedDeltaTracking; 	

		if (std::string(argv[i]) == "-ratiotracking")
			method = method | TrackingMethod::TMRatioTracking; 	

		if (std::string(argv[i]) == "-residualratiotracking")
			method = method | TrackingMethod::TMResidualRatioTracking; 	

		if (std::string(argv[i]) == "-quadrature")
			method = method | TrackingMethod::TMQuadrature; 	

		if (std::string(argv[i]) == "-all")
			method = method | TrackingMethod::TMAll; 	

	}

	if(!method) method = method | TrackingMethod::TMRayMarching;

	auto [func,gt] = function1d(argc,argv);
	
	svg_cpp_plot::SVG svg;
	svg.viewBox(svg_cpp_plot::BoundingBox(-60,-60,1720,860));


	// Plot Function
	svg_cpp_plot::_2d::polyline plot_function;

	double f_min_val = 1., f_max_val = -1., f_exp_value = gt(0.f,1.f);

	for(double s=0.; s<1.; s+=1.e-4)
	{
		auto fs = func(s);
		plot_function.add_point(s,fs);

		if(fs < f_min_val) f_min_val = fs;
		if(fs > f_max_val) f_max_val = fs;

	}

	auto& graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(	svg_cpp_plot::Graph2D({800,800},
																							svg_cpp_plot::BoundingBox(0,f_min_val,double(1),f_max_val)));
	graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xticks(10).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.yticks(10).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xlabels(4);
	graph.ylabels(10);

	graph.area().add(plot_function).stroke(svg_cpp_plot::red).stroke_width(2);

	svg_cpp_plot::_2d::polyline plot_rm, plot_dt, plot_wdt, plot_rt, plot_rrt, plot_quad;




	printf("Ground Truth: %.10f - Max Transmittance:  %.10f -- Min Transmittance: %.1f\n", exp(-f_exp_value), f_max_val, f_min_val);


	double min_val = 1., max_val = -1.;


	for(unsigned int i = 1; i<max_log2_samples; ++i)
	{

		unsigned long sp = pow(2,i);



		if( method & TrackingMethod::TMRayMarching )
		{
			// Ray Marching
			auto v_raymarching(ray_marching(sp));
			using ResultRM = decltype(v_raymarching.init(func, range(0.f,1.f)));
			ResultRM v_rmsamples;

			//const F& f, const Range<Float,1>& range, Samples<Result>& samples
			v_raymarching.step(func, range(0.f,1.f), v_rmsamples, verbose);
			float rmintegral = v_raymarching.integral(func, v_rmsamples);
			unsigned long rmnb_queries = v_raymarching.queries(func, v_rmsamples);

			plot_rm.add_point(i,rmintegral);

			if(rmintegral < min_val && rmintegral > 0) min_val = rmintegral;
			if(rmintegral > max_val) max_val = rmintegral;

			printf("RayMarching\t- Integral: %.10f [%d queries] \n", rmintegral, rmnb_queries);
		}

		// Delta Tracking
		if( method & TrackingMethod::TMDeltaTracking )
		{
			auto v_deltatracking(delta_tracking(std::mt19937_64(seed),f_max_val));

			using ResultDeltaTracking = decltype(v_deltatracking.init(func, range(0.f,1.f)));
			ResultDeltaTracking v_dtsamples;

			for(unsigned long j=0; j<sp; ++j)
				v_deltatracking.step(func, range(0.f,1.f), v_dtsamples, verbose);
		
			float dtintegral = v_deltatracking.integral(func, v_dtsamples);
			unsigned long dtnb_queries = v_deltatracking.queries(func, v_dtsamples);

			if(dtintegral < min_val && dtintegral > 0) min_val = dtintegral;
			if(dtintegral > max_val) max_val = dtintegral;
			plot_dt.add_point(i,dtintegral);
			printf("DeltaTracking\t- Integral: %.10f [%d samples - %d queries] \n", dtintegral, sp, dtnb_queries);
		}	

		// Weighted Delta Tracking
		if( method & TrackingMethod::TMWeightedDeltaTracking )
		{
			auto v_weighteddeltatracking(weighted_delta_tracking(std::mt19937_64(seed),f_exp_value));

			using ResultType = decltype(v_weighteddeltatracking.init(func, range(0.f,1.f)));
			ResultType v_samples;

			for(unsigned long j=0; j<sp; ++j)
				v_weighteddeltatracking.step(func, range(0.f,1.f), v_samples, verbose);
		
			float integral = v_weighteddeltatracking.integral(func, v_samples);
			unsigned long nb_queries = v_weighteddeltatracking.queries(func, v_samples);

			if(integral < min_val && integral > 0) min_val = integral;
			if(integral > max_val) max_val = integral;
			plot_wdt.add_point(i,integral);
			printf("WgtDeltaTracking\t- Integral: %.10f [%d samples - %d queries] \n", integral, sp, nb_queries);
		}	


		// Ratio Tracking
		if( method & TrackingMethod::TMRatioTracking )
		{
			auto v_ratiotracking(ratio_tracking(std::mt19937_64(seed),f_max_val));
			using ResultRatTracking = decltype(v_ratiotracking.init(func, range(0.f,1.f)));
			ResultRatTracking v_rtsamples;

			for(unsigned long j=0; j<sp; ++j)
				v_ratiotracking.step(func, range(0.f,1.f), v_rtsamples, verbose);
		
			float rtintegral = v_ratiotracking.integral(func, v_rtsamples);
			unsigned long rtnb_queries = v_ratiotracking.queries(func, v_rtsamples);


			if(rtintegral < min_val && rtintegral > 0) min_val = rtintegral;
			if(rtintegral > max_val) max_val = rtintegral;
			plot_rt.add_point(i,rtintegral);
			printf("RatioTracking\t- Integral: %.10f [%d samples - %d queries] \n", rtintegral, sp, rtnb_queries);
		}

		// Residual Ratio Tracking
		if( method & TrackingMethod::TMResidualRatioTracking )
		{
			auto v_residualratiotracking(residual_ratio_tracking(std::mt19937_64(seed),f_max_val, f_exp_value));
			using ResultResRatTracking = decltype(v_residualratiotracking.init(func, range(0.0,1.0)));
			ResultResRatTracking v_rrtsamples;

			for(unsigned long j=0; j<sp; ++j)
				v_residualratiotracking.step(func, range(0.,1.), v_rrtsamples, verbose);
		
			float rrtintegral = v_residualratiotracking.integral(func, v_rrtsamples);
			unsigned long rrtnb_queries = v_residualratiotracking.queries(func, v_rrtsamples);


			if(rrtintegral < min_val && rrtintegral > 0) min_val = rrtintegral;
			if(rrtintegral > max_val) max_val = rrtintegral;
			plot_rrt.add_point(i,rrtintegral);
			printf("RRatioTracking\t- Integral: %.10f [%d samples - %d queries] \n", rrtintegral, sp, rrtnb_queries);
		}

		// Quadrature
		if( method & TrackingMethod::TMQuadrature )
		{
			auto v_method(stepper_adaptive(nested(simpson,trapezoidal)));
			FunctionWrapper<decltype(func),1> v_F(func);


			using StepperData = decltype(v_method.init(v_F, range(0.,1.)));
			std::vector<StepperData> v_samples;


			v_samples.push_back(v_method.init(v_F, range(0.,1.)));

			for(unsigned long j=0; j<sp; ++j)
				//for (auto& d : v_samples) 
				v_method.step(v_F, range(0.,1.), v_samples[0]);
		
			double integral = exp(-v_method.integral(v_F, range(0.,1.), v_samples[0]));
			unsigned long nb_queries = v_F.samples();


			if(integral < min_val && integral > 0) min_val = integral;
			if(integral > max_val) max_val = integral;
			plot_quad.add_point(i,integral);
			printf("SimpsonTrapez\t- Integral: %.10f [%d samples - %d queries] \n", integral, sp, nb_queries);
		}
	}

	printf("Ground Truth\t- %.10f \n", exp(-gt(0.f,1.f)));
//	printf("Min Val: %f - Max Val: %f\n", min_val, max_val);

	auto graph1 = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(860,0))).add(	svg_cpp_plot::Graph2D({800,800},
																							svg_cpp_plot::BoundingBox(0,0,max_log2_samples,1)));

	graph1.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph1.xticks(10).stroke(svg_cpp_plot::black).stroke_width(1);
	graph1.yticks(10).stroke(svg_cpp_plot::black).stroke_width(1);
	graph1.xlabels(4);
	graph1.ylabels(10);

	graph1.area().add(plot_rm).stroke(svg_cpp_plot::red).stroke_width(2);
	graph1.area().add(plot_dt).stroke(svg_cpp_plot::green).stroke_width(2);
	graph1.area().add(plot_wdt).stroke(svg_cpp_plot::green).stroke_width(4);
	graph1.area().add(plot_rt).stroke(svg_cpp_plot::blue).stroke_width(2);
	graph1.area().add(plot_rrt).stroke(svg_cpp_plot::blue).stroke_width(4);
	graph1.area().add(plot_quad).stroke(svg_cpp_plot::orange).stroke_width(2);

	std::ofstream saving_file(output);
	saving_file << svg;


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



