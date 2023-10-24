#include "../functions/functions-render3d.h"
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
	
	auto render_function = std::get<0>(render_function3d(argc,argv));
	std::array<float,3> range_min, range_max; range_min.fill(0); range_max.fill(1);
	Range<float,3> range(range_min, range_max);
	
    CImgWrapper<float> image(w,h);
    tracer::Scene scene;
   
    image.clear();
    integrate_bins_stepper_progression("Monte-Carlo             ",stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),spp,image,image.resolution(),render_function, range);
    image.save(output+"-monte-carlo.hdr");
    image.clear();
    integrate_bins_stepper_progression("Simpson-Trapezoidal     ",stepper_bins_adaptive(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5)),(spp*w*h)/(3*3*2),image,image.resolution(),render_function, range);
    image.save(output+"-simpson-trapezoidal.hdr");
    image.clear();
    integrate_bins_stepper_progression("Simpson-Trapezoidal CV  ",stepper_bins_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(1.e-5), (spp*w*h)/(3*3*4), seed, seed+1),spp/2,image,image.resolution(),render_function, range);
    image.save(output+"-simpson-trapezoidal-cv.hdr");
}




