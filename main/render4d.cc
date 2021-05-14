#include "../functions/functions-render4d.h"
#include "../quadrature/range.h"
#include "../quadrature/integrate-adaptive-control-variates.h"
#include "../quadrature/integrate-bins-stepper.h"
#include "../quadrature/integrate-bins-adaptive.h"
#include "../quadrature/integrate-bins-adaptive-precalculate.h"
#include "../utils/cimg-wrapper.h"
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>

int main(int argc, char** argv) {
    int w = 512; int h = 512;
	unsigned long spp = 128;
    std::size_t seed = std::random_device()();
	std::string output = "output";
	
	for (int i = 0; i<(argc-1); ++i) {
		if (std::string(argv[i])=="-width") w = atoi(argv[++i]);
		else if (std::string(argv[i])=="-height") h = atoi(argv[++i]);
		else if (std::string(argv[i])=="-spp") spp = atol(argv[++i]);
        else if (std::string(argv[i])=="-seed") seed = atol(argv[++i]);
		else if (std::string(argv[i])=="-output") output = std::string(argv[++i]);
	}
	
	auto render_function = std::get<0>(render_function4d(argc,argv));
	std::array<float,4> range_min, range_max; range_min.fill(0); range_max.fill(1);
	Range<float,4> range(range_min, range_max);
	
    CImgWrapper<float> image(w,h);
    tracer::Scene scene;
   
    image.clear();
    integrate_bins_stepper_progression("Monte-Carlo             ",stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),spp,image,image.resolution(),render_function, range);
    image.save(output+"-monte-carlo.hdr");
    image.clear();
    integrate_bins_stepper_progression("Simpson-Trapezoidal     ",stepper_bins_adaptive(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5)),(spp*w*h)/(3*3*3*2),image,image.resolution(),render_function, range);
    image.save(output+"-simpson-trapezoidal.hdr");
	/*
	image.clear();
	integrate_bins_stepper("Prec.Simpson-Trapezoidal",stepper_bins_adaptive_precalculate(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5)),(spp*w*h)/(3*3*3*2),image,image.resolution(),render_function, render_function.range());
    image.save("render-area-simpson-trapezoidal-precalculate.hdr");
	*/
    image.clear();
    integrate_bins_stepper_progression("Boole-Simpson           ",stepper_bins_adaptive(nested(boole,simpson)),(spp*w*h)/(5*5*5*4),image,image.resolution(),render_function, range);
    image.save(output+"-boole-simpson.hdr");

	/*
	image.clear();
    integrate_bins_stepper("Prec.Boole-Simpson      ",stepper_bins_adaptive_precalculate(nested(boole,simpson)),(spp*w*h)/(5*5*5*4),image,image.resolution(),render_function, render_function.range());
    image.save("render-area-boole-simpson-precalculate.hdr");
	*/
    image.clear();
    integrate_bins_stepper_progression("Simpson-Trapezoidal CV  ",stepper_bins_adaptive_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5),stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),vector_sampler_uniform(seed+1),(spp*w*h)/(3*3*3*4)),(spp*w*h)/(3*3*3*4),image,image.resolution(),render_function, range);
    image.save(output+"-simpson-trapezoidal-cv.hdr");
}




