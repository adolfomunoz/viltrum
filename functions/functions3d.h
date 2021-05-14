#pragma once

#include <functional>
#include "../quadrature/monte-carlo.h"
#include <PerlinNoise/PerlinNoise.hpp>
#include <cimg-all.h>

using Function3D = std::tuple<std::function<double(double,double,double)>,std::function<double(double,double,double,double,double,double)>>;

Function3D function3d(int& i, int argc, char **argv) {
    std::function<double(double,double,double)> func = [] (double x, double y, double z) { return 0.5f; };
	std::function<double(double,double,double,double,double,double)> groundtruth;
    unsigned long ground_truth_samples = 128000;
    std::size_t seed = std::random_device()();
    bool has_groundtruth = false;
    for (int j = 0; j < argc-1; ++j) {
        if (std::string(argv[j])=="-ground-truth-samples") ground_truth_samples = atol(argv[++j]);
        else if (std::string(argv[j])=="-ground-truth-seed") seed = atol(argv[++j]);
    }
    
    if (std::string(argv[i])=="cartesian2x1") {
        auto f1 = std::get<0>(function2d(++i, argc, argv));
        auto f2 = std::get<0>(function1d(++i, argc, argv));
        func = [f1,f2] (double x, double y, double z) {     return f1(x,y)*f2(z); };
    }
    else if (std::string(argv[i])=="interpolation") {
        auto f1 = std::get<0>(function2d(++i, argc, argv));
        auto f2 = std::get<0>(function2d(++i, argc, argv));
        func = [f1,f2] (double x, double y, double z) {     return f1(x,y)*z + (1-z)*f2(x,y); };
        //We could calculate the ground truth here based
        //on f1 and f2 ground truth. Maybe will do it
        //eventually
    }
    else if (std::string(argv[i]) == "perlin") {
        double freq(1); if ((i<(argc-1)) && (argv[i+1][0]!='-')) freq=atof(argv[++i]);
        std::size_t seed = std::random_device()(); 
        if ((i<(argc-1)) && (argv[i+1][0]!='-')) seed=atol(argv[++i]);
        func = [freq,seed] (double x, double y, double z) { return 0.5+0.5*siv::PerlinNoise(std::uint32_t(seed)).noise3D(x*freq,y*freq,z*freq); };
    }
    else if (std::string(argv[i]) == "video") {
        // load_ffmpeg_external
        // get_load_ffmpeg_external
        cimg_library::CImg<double> input(argv[++i]);
        cimg_library::CImg<double> grayscale = input.get_channel(0);
        for (int c=1;c<input.spectrum();++c) grayscale += input.get_channel(c);
        grayscale /= grayscale.max();
        func = [grayscale] (double x, double y, double z)
             { return grayscale.cubic_atXYZ(
                float(x)*float(grayscale.width()),
                float(1.0f-y)*float(grayscale.height()),std::max(0.0f,std::min(float(grayscale.depth()-1),float(z)*float(grayscale.depth()))),0); };
    }

    if (!has_groundtruth) groundtruth = [func,ground_truth_samples,seed] (double a0, double a1, double a2, double b0, double b1, double b2) {
        return integrator_monte_carlo_uniform(ground_truth_samples,seed).integrate([func] (const std::array<double,3>& x) { return func(std::get<0>(x),std::get<1>(x), std::get<2>(x)); },range(a0,a1,a2,b0,b1,b2));
	};

	return std::make_tuple(func,groundtruth);     
}

Function3D function3d(int argc, char **argv) {
	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-function") 
            return function3d(++i,argc,argv);
	}
    
    return Function3D(
        [] (double x, double y, double z) { return 0.5f; },
        [] (double ax, double ay, double az, double bx, double by, double bz) { return 0.5f*(bx-ax)*(by-ay)*(bz-az);});
}

