#pragma once

#include <functional>
#include "../quadrature/monte-carlo.h"
#include <PerlinNoise/PerlinNoise.hpp>
#include "BrownianNoise.hpp"

using Function1D = std::tuple<std::function<double(double)>,std::function<double(double,double)>>;

Function1D function1d(int& i, int argc, char** argv) {
    std::function<double(double)> f = [] (double x) { return 0.5f; };
	std::function<double(double,double)> groundtruth;
    unsigned long ground_truth_samples = 10000000;
    std::size_t seed = std::random_device()();
    bool has_groundtruth = false;
    for (int j = 0; j < argc-1; ++j) {
        if (std::string(argv[j])=="-ground-truth-samples") ground_truth_samples = atol(argv[++j]);
        else if (std::string(argv[j])=="-ground-truth-seed") seed = atol(argv[++j]);
    }
    
    std::string function = argv[i];

    if (function == "cos") {
        double param = atof(argv[++i]);
        f = [param] (double x) { return std::cos(2*M_PI*x/param); };
        groundtruth = [param] (double a, double b) {
            return param*(std::sin(2*M_PI*b/param) - std::sin(2*M_PI*a/param))/(2*M_PI);
        };
        has_groundtruth=true;
    }
    if (function == "abscos") {
        double param = atof(argv[++i]);
        f = [param] (double x) { return std::cos(2*M_PI*x/param)+1; };
        groundtruth = [param] (double a, double b) {
            return b-a + param*(std::sin(2*M_PI*b/param) - std::sin(2*M_PI*a/param))/(2*M_PI);
        };
        has_groundtruth=true;
    }
    else if (function == "step") {
        double param(0); if ((i<(argc-2)) && std::string(argv[i+1]) == "-param"){ param=atof(argv[i+2]); i+=2;}
        double offset(0); if ((i<(argc-2)) && std::string(argv[i+1]) == "-offset"){ offset=atof(argv[i+2]); i+=2;}

        f = [param, offset] (double x) { return ((x<param)?-1.0:1.0)+offset; };
        groundtruth = [param, offset] (double a, double b) {
            return (b-a)*offset + 	(((a<param)&&(b<param))?(a-b):
                                    ((a>param)&&(b>param))?(b-a):
                                                           (a+b - 2*param));
        };

        has_groundtruth=true;
    }
    else if (function == "id") {
        f = [] (double x) { return x; };
    }
    else if (function == "perlin") {
        double freq(1); if ((i<(argc-2)) && std::string(argv[i+1]) == "-frequency"){ freq=atof(argv[i+2]); i+=2;}
        double offset(0); if ((i<(argc-2)) && std::string(argv[i+1]) == "-offset"){ offset=atof(argv[i+2]); i+=2;}
        std::size_t seed = std::random_device()(); 
        if ((i<(argc-1)) && (argv[i+1][0]!='-')) seed=atol(argv[++i]);

        f = [freq,seed, offset] (double x) { return siv::PerlinNoise(std::uint32_t(seed)).noise1D(x*freq) + offset; };
    }
    else if (function == "fbm") {
        double freq(1); if ((i<(argc-2)) && std::string(argv[i+1]) == "-frequency"){ freq=atof(argv[i+2]); i+=2;}
        double offset(0); if ((i<(argc-2)) && std::string(argv[i+1]) == "-offset"){ offset=atof(argv[i+2]); i+=2;}
        std::size_t seed = std::random_device()(); 
        f = [freq,seed,offset] (double x) { return Noise::BrownianNoise().noise(x*freq) + offset; };
    }   

    if (!has_groundtruth) groundtruth = [f,ground_truth_samples,seed] (double a, double b) {
        return integrator_monte_carlo_uniform(ground_truth_samples,seed).integrate([f] (const std::array<double,1>& x) { return f(std::get<0>(x)); },range(a,b));
	};
	return std::make_tuple(f,groundtruth);    
}

std::tuple<std::function<double(double)>,std::function<double(double,double)>> function1d(int argc, char **argv) {
	for (int i = 0; i < argc; ++i) {
		if (std::string(argv[i]) == "-function")   
            return function1d(++i,argc,argv);
	}
    return std::tuple<std::function<double(double)>,std::function<double(double,double)>>(
        [] (double x) { return 0.5f; },
        [] (double a, double b) { return 0.5f*(b-a);});
}

