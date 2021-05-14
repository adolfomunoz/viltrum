#pragma once

#include <functional>
#include "../quadrature/monte-carlo.h"
#include <PerlinNoise/PerlinNoise.hpp>
#include <cimg-all.h>

using Function2D = std::tuple<std::function<double(double,double)>,std::function<double(double,double,double,double)>>;

Function2D function2d(int& i, int argc, char **argv) {
    std::function<double(double,double)> func = [] (double x,double y) { return 0.5f; };
	std::function<double(double,double,double,double)> groundtruth;
    unsigned long ground_truth_samples = 128000;
    std::size_t seed = std::random_device()();
    bool has_groundtruth = false;
    for (int j = 0; j < argc-1; ++j) {
        if (std::string(argv[j])=="-ground-truth-samples") ground_truth_samples = atol(argv[++j]);
        else if (std::string(argv[j])=="-ground-truth-seed") seed = atol(argv[++j]);
    }
    
    if (std::string(argv[i])=="step") {
        double atx = atof(argv[++i]);
        double aty = atof(argv[++i]);
        func = [atx,aty] (double x, double y) { return (((x/atx)+(y/aty))<1.0)?0.1:0.9; };
    }
    else if (std::string(argv[i])=="sin") {
        double freqx = atof(argv[++i]);
        double freqy = atof(argv[++i]);
        func = [freqx,freqy] (double x, double y) { return 0.5+0.4*std::sin(2*M_PI*x*freqx)*std::sin(2*M_PI*y*freqy); };
    }
    else if (std::string(argv[i])=="freqband") {
        double freq = atof(argv[++i]);
        func = [freq] (double x, double y) {
            return 0.5*y + (1-y)*(0.5+0.5*std::sin(2*M_PI*std::exp(x*freq)));
        };
    }
    else if (std::string(argv[i])=="lightmedia") {
        double extinction = atof(argv[++i]);
        func = [extinction] (double x, double y) {
            double d = std::sqrt((x-0.5)*(x-0.5) + (1.01-y)*(1.01-y));
            return std::exp(-extinction*(x+d))/(d*d);
        };
    }					
    else if (std::string(argv[i])=="box") {
        double at1 = atof(argv[++i]); 
        double at2 = atof(argv[++i]);
        func = [at1,at2] (double x, double y) { return (x<at1)?0.2:((y>at2)?0.7:0.8); };
    }			
    else if (std::string(argv[i])=="horizontal") {
        double freq = atof(argv[++i]);
        func = [freq] (double x, double y) { return 0.5+0.4*std::sin(2*M_PI*x*freq); };
    }
    else if (std::string(argv[i])=="step-smooth") {
        double at = atof(argv[++i]);
        func = [at] (double x, double y) { return (x<at)?(0.1+y*0.4):(0.9-y*0.4); };
    }
    else if (std::string(argv[i]) == "perlin") {
        double freq(1); if ((i<(argc-1)) && (argv[i+1][0]!='-')) freq=atof(argv[++i]);
        std::size_t seed = std::random_device()(); 
        if ((i<(argc-1)) && (argv[i+1][0]!='-')) seed=atol(argv[++i]);
        func = [freq,seed] (double x, double y) { return 0.5+0.5*siv::PerlinNoise(std::uint32_t(seed)).noise2D(x*freq,y*freq); };
    }
    else if (std::string(argv[i]) == "image") {
        cimg_library::CImg<double> input(argv[++i]);
        cimg_library::CImg<double> grayscale = input.get_channel(0);
        for (int c=1;c<input.spectrum();++c) grayscale += input.get_channel(c);
        grayscale /= grayscale.max();
        func = [grayscale] (double x, double y)
             { return grayscale.cubic_atXY(float(x)*float(grayscale.width()),float(1.0f-y)*float(grayscale.height()),0,0); };
    }

    if (!has_groundtruth) groundtruth = [func,ground_truth_samples,seed] (double a0, double a1, double b0, double b1) {
        return integrator_monte_carlo_uniform(ground_truth_samples,seed).integrate([func] (const std::array<double,2>& x) { return func(std::get<0>(x),std::get<1>(x)); },range(a0,a1,b0,b1));
	};

	return std::make_tuple(func,groundtruth);     
}

std::tuple<std::function<double(double,double)>,std::function<double(double,double,double,double)>> function2d(int argc, char **argv) {
	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-function") 
            return function2d(++i,argc,argv);
	}
    
    return Function2D(
        [] (double x, double y) { return 0.5f; },
        [] (double ax, double ay, double bx, double by) { return 0.5f*(bx-ax)*(by-ay);});
}

