#include <memory>
#include <svg-cpp-plot/svg-cpp-plot.h>
#include "../quadrature/integrate-bins.h"
#include "../quadrature/integrate-bins-stepper.h"
#include "../quadrature/monte-carlo.h"
#include "../quadrature/control-variates.h"
#include "../plot/convergence.h"
#include "../plot/control-variate.h"
#include "../functions/functions2d.h"
#include "../quadrature/integrate-adaptive-control-variates.h"
#include "../quadrature/integrate-optimized-adaptive-stratified-control-variates.h"

#include <iostream>
#include <cmath>

int main(int argc, char **argv) {	
	const char* output = "output.svg";
	double max_samples_or_time = 1000;
	bool plot_samples = true;
	std::size_t seed = std::random_device()();
	auto [func,gt] = function2d(argc,argv);
//	double rmin = 0; double rmax = 1;
	double min_plottable_error = 1.e-10;
	int plot_resolution = 100;
    int plot_image_resolution = 300;
    std::size_t montecarlo_error_average = 100;
    unsigned long adaptive_iterations = 11;
	float cv_ratio = 1.0f/16.0f;
    unsigned int ncontrolvariates = 0;
	std::size_t bins = 9;
    double error_size_weight = 1.e-3;

    bool test_order = false;
    bool test_errorstrategy = false;
    bool relative_cv = false;
	bool test_quadrature = false;

    int width = 400;
    int spacing = 40;


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
		else if (std::string(argv[i])=="-plot-image-resolution") { plot_image_resolution = atoi(argv[++i]); }
		else if (std::string(argv[i])=="-adaptive-iterations") { if (ncontrolvariates==0) ncontrolvariates=1; adaptive_iterations = atol(argv[++i]); }
		else if (std::string(argv[i])=="-cv-ratio") { if (ncontrolvariates==0) ncontrolvariates=1; cv_ratio=atof(argv[++i]); relative_cv=true; }
        else if (std::string(argv[i])=="-number-of-cv-tests") { ncontrolvariates=atoi(argv[++i]); } 
        else if (std::string(argv[i])=="-error-size-weight") { error_size_weight=atof(argv[++i]); }
		else if (std::string(argv[i])=="-bins") { bins = atoi(argv[++i]); }         
		else if (std::string(argv[i])=="-width") { width = atoi(argv[++i]); } 
		else if (std::string(argv[i])=="-spacing") { spacing = atoi(argv[++i]); } 
	}

	for (int i = 0; i<argc; ++i) {
		if (std::string(argv[i])=="-test-order") { test_order = true; } 
        else if (std::string(argv[i])=="-test-errorstrategy") { test_errorstrategy = true; } 
        else if (std::string(argv[i])=="-integrator-modern") { relative_cv=true; }
        else if (std::string(argv[i])=="-relative-cv") { relative_cv=true; }
        else if (std::string(argv[i])=="-test-quadrature") { test_quadrature=true; }

	}


//    auto color_mc  = svg_cpp_plot::rgb(43.0/255.0,109.0/255.0,167.0/255.0);
    auto color_mc  = svg_cpp_plot::rgb(1,0,0);
    auto color_ours  = svg_cpp_plot::rgb(30.0/255.0,149.0/255.0,93.0/255.0);
    float graph_width = 2;

    std::vector<double> ground_truth(bins);
    double dx = 1.0/double(bins);
    for (std::size_t i=0; i<bins; ++i) ground_truth[i]=double(bins)*gt(dx*i,0,dx*(i+1),1);
	


    auto error = [ground_truth,min_plottable_error] (const std::vector<double>& a) {
        double err = 0;
        for (std::size_t i=0;i<a.size();++i)
            err += std::abs(ground_truth[i]-a[i]);
		return std::max(min_plottable_error,err/double(a.size()));
	};

    std::cerr<<"Sampling error average = "<<montecarlo_error_average<<std::endl;

	svg_cpp_plot::SVG svg;
    	svg.viewBox(svg_cpp_plot::BoundingBox(-2*spacing,-spacing,width+spacing,width+2*spacing));

    std::list<svg_cpp_plot::_2d::polyline> plots;
	std::cout<<"Monte-Carlo"<<std::endl;
    plots.push_back(convergence_bins(plot_samples,stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, bins, montecarlo_error_average, plot_resolution)
		.stroke(color_mc).stroke_width(graph_width));
    auto [min_x, max_error] = plots.back().point_list().front(); // We set the maximum error and the minimum x to the maximum Monte-Carlo error

    if (test_quadrature) {
		std::cout<<"Simpson-Trapezoidal"<<std::endl;  
		plots.push_back(convergence_bins(plot_samples,stepper_bins_adaptive(nested(simpson,trapezoidal),error_single_dimension_size(error_size_weight)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time,bins,1,plot_resolution)
		    .stroke(svg_cpp_plot::rgb(0,1,0)).stroke_width(graph_width));
    }

    if (plot_samples && relative_cv && (ncontrolvariates>0)) {
		float ratio; unsigned int i;
		for (i=0, ratio=cv_ratio; i<ncontrolvariates; ++i, ratio*=cv_ratio) {
            float dc = 2.0/float(ncontrolvariates+2);        
			std::cout<<"Control Variate at 1/"<<(1.0/ratio)<<std::endl;
			plots.push_back(convergence_bins_integrator_samples([&] (unsigned long spp) {
                    unsigned long spp_cv = std::max(1UL,(unsigned long)(ratio*spp));
                    return integrator_optimized_perregion_adaptive_stratified_control_variates(nested(simpson, trapezoidal),error_single_dimension_size(error_size_weight),(unsigned long)(bins*spp_cv/(3*2)),std::max(1L,long(spp)-long(spp_cv)),seed);
                    },func, range(0.0,0.0,1.0,1.0),error,max_samples_or_time, bins, montecarlo_error_average, plot_resolution).stroke(svg_cpp_plot::rgb(std::min(2.0f-(i+1)*dc,1.0f),std::min((i+1)*dc,1.0f),0)).stroke_width(graph_width));
		}
    }


    if (test_order && test_quadrature) {
		std::cout<<"Boole-Simpson"<<std::endl;  
        plots.push_back(convergence_bins(plot_samples,stepper_bins_adaptive(nested(boole,simpson),error_single_dimension_size(error_size_weight)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time,bins,1,plot_resolution)
        	.stroke(svg_cpp_plot::rgb(0,0,1)).stroke_width(2));
	}
    
    if (test_errorstrategy && test_quadrature) {
		std::cout<<"Simpson-Trapezoidal + error strategy"<<std::endl;  
        plots.push_back(convergence_bins(plot_samples,stepper_bins_adaptive(nested(simpson,trapezoidal)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time,bins,1,plot_resolution)
                .stroke(svg_cpp_plot::rgb(0,1,0)).stroke_width(2).stroke_dasharray({2,4}));
	}

    if (test_errorstrategy && test_order && test_quadrature) {
		std::cout<<"Boole-Simpson + error strategy"<<std::endl;  
        plots.push_back(convergence_bins(plot_samples,stepper_bins_adaptive(nested(boole,simpson)), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time,bins,1,plot_resolution)
		    .stroke(svg_cpp_plot::rgb(0,0,1)).stroke_width(graph_width).stroke_dasharray({graph_width,2*graph_width}));
	}
	
    if ((ncontrolvariates>0) && (!relative_cv)) {
        for (unsigned int i = 1; i<=ncontrolvariates; ++i) {
            unsigned long iterations = std::pow(adaptive_iterations,i);
		    std::cout<<"Simpson-Trapezoidal + control variate "<<iterations<<" iterations"<<std::endl;
            float dc = 2.0/float(ncontrolvariates+2);            
			plots.push_back(convergence_bins_integrator_samples([&] (unsigned long spp) {
					unsigned long spp_residual = float(bins*spp - iterations*3*2)/float(bins); 
                    return integrator_optimized_perregion_adaptive_stratified_control_variates(nested(simpson, trapezoidal),error_single_dimension_size(error_size_weight),iterations,std::max(1UL,spp_residual),seed);
                    },func, range(0.0,0.0,1.0,1.0),error,max_samples_or_time, bins, montecarlo_error_average, plot_resolution).stroke(svg_cpp_plot::rgb(std::min(2.0f-(i+1)*dc,1.0f),std::min((i+1)*dc,1.0f),0)).stroke_width(graph_width));
		}
	}

    if ((ncontrolvariates>0) && test_order&& (!relative_cv)) {
        for (unsigned int i = 1; i<=ncontrolvariates; ++i) {
            unsigned long iterations = std::pow(adaptive_iterations,i);
		    std::cout<<"Boole-Simpson + control variate "<<iterations<<std::endl;  		
            float dc = 2.0/float(ncontrolvariates+1);        
            plots.push_back(convergence_bins(plot_samples,stepper_bins_adaptive_stratified_control_variates(nested(boole, simpson),error_single_dimension_size(error_size_weight),iterations, seed, seed+1), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average,plot_resolution)
		        .stroke(svg_cpp_plot::rgb(std::min(2.0f-i*dc,1.0f),0,std::min(i*dc,1.0f))).stroke_width(2));
        }
	}

    float min_error = max_error; float max_x = std::get<0>(plots.front().point_list().back());
    for (const auto& p: plots) {
        auto [x, err] = p.point_list().back(); // We check the last point of each plot
        if (err<min_error) min_error=err;
        if (x < max_x) max_x = x;
    }

	auto& graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(0,0))).add(svg_cpp_plot::Graph2D({width,width},svg_cpp_plot::BoundingBox(min_x,min_error,max_x,max_error)));
	graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
	graph.xticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
	graph.yticks(2).stroke(svg_cpp_plot::black).stroke_width(1);
    for (const auto& p: plots) graph.area().add(p);
//	std::string xlabel = plot_samples?"log(samples)":"log(time)";
//	graph.add(svg_cpp_plot::_2d::text({width/2,width+14},xlabel)).font_size(30).text_anchor(svg_cpp_plot::text_anchor_middle).dominant_baseline(svg_cpp_plot::dominant_baseline_hanging);
//	using svg_cpp_plot::_2d::operator*;

//	graph.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(-14,width/2)*svg_cpp_plot::_2d::rotate(-M_PI/2))).add(svg_cpp_plot::_2d::text({0,0},"log(error)")).font_size(30).text_anchor(svg_cpp_plot::text_anchor_middle);
//	graph.xlabels(2);
//	graph.ylabels(2);

	/*For the function plot to be inside is should be in (540,10) and (540,260) but we are getting it out*/
/*	auto& function_graph = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(820,0))).add(svg_cpp_plot::Graph2D({250,250},svg_cpp_plot::BoundingBox(0,0,1,1)));
	function_graph.border().stroke(svg_cpp_plot::black).stroke_width(1);
    function_graph.area().add(svg_cpp_plot::_2d::function_2d(func,{0.0,0.0},{1.0,1.0},{plot_image_resolution,plot_image_resolution}));
	auto& function_gt = svg.add(svg_cpp_plot::_2d::group(svg_cpp_plot::_2d::translate(820,250))).add(svg_cpp_plot::Graph2D({250,100},svg_cpp_plot::BoundingBox(0,0,1,1)));
	function_gt.border().stroke(svg_cpp_plot::black).stroke_width(1);
    function_gt.area().add(svg_cpp_plot::_2d::bar_plot(ground_truth)).stroke_width(0.5).stroke(svg_cpp_plot::black).fill(svg_cpp_plot::yellow);
*/	
	std::ofstream f(output);
	f << svg;
	/*
    graph.area().add(convergence_bins(plot_samples,stepper_bins_control_variate(control_variate_quadrature_adaptive(nested(simpson, trapezoidal),adaptive_iterations), seed), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0,1,0)).stroke_width(2);
    graph.area().add(convergence_bins(plot_samples,stepper_bins_control_variate(control_variate_quadrature_adaptive(nested(boole, simpson),adaptive_iterations), seed), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0,0,1)).stroke_width(2);
    graph.area().add(convergence_bins(plot_samples,stepper_bins_adaptive_control_variates(nested(simpson, trapezoidal),adaptive_iterations, seed, seed+1), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0.5,1,0.5)).stroke_width(2);
    graph.area().add(convergence_bins(plot_samples,stepper_bins_adaptive_control_variates(nested(boole, simpson),adaptive_iterations, seed, seed+1), func, range(0.0,0.0,1.0,1.0), error, max_samples_or_time, montecarlo_error_average))
		.stroke(svg_cpp_plot::rgb(0.5,0.5,1)).stroke_width(2);
*/

	return 0;
}



