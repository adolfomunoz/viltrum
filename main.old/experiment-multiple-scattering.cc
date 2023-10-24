#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "../render/multiple-scattering.h"
#include "../quadrature/range.h"
#include "../quadrature/integrate-adaptive-control-variates.h"
#include "../quadrature/integrate-bins-stepper.h"
#include "../quadrature/integrate-bins-adaptive.h"
#include "../quadrature/integrate-bins-adaptive-precalculate.h"
#include "../utils/cimg-wrapper.h"



Medium medium_from_commandline(int argc, char** argv) {
    float absorption = 0.1f, scattering = 0.9f;
    for (int i = 0; i<(argc-1);++i) {
        if (std::string(argv[i])=="-absorption") absorption = atof(argv[++i]);
        else if (std::string(argv[i])=="-scattering") scattering = atof(argv[++i]);
    }
    return Medium(Spectrum::Constant(absorption),Spectrum::Constant(scattering));
}

tracer::Scene scene_from_commandline(int argc, char** argv) {
    tracer::Scene scene; 
    float scale = 1.0f;
    for (int i = 0; i<(argc-1);++i) {
        if (std::string(argv[i])=="-scale") scale = atof(argv[++i]);
    }

    for (int i = 0;i<argc;++i) {
        if (std::string(argv[i])=="occluder")   scene.add(tracer::Sphere(Eigen::Vector3f(0.5*scale,0,0),0.2*scale).set_material(lambertian(Spectrum::Constant(0.7))));
        else if (std::string(argv[i])=="plane") scene.add(tracer::Plane(Eigen::Vector3f(1,0,0),Eigen::Vector3f(-1*scale,0,0)).set_material(lambertian(Spectrum::Constant(0.7))));
		else if (std::string(argv[i])=="arealight")   scene.add(tracer::Sphere(Eigen::Vector3f(1.0*scale,0,0),0.1*scale).set_material(emitter(Spectrum::Constant(100.0))));
    }
    return scene;
}

tracer::Pinhole camera_from_commandline(int argc, char** argv) {
    float scale = 1.0f; int w = 512; int h = 512;
    for (int i = 0; i<(argc-1);++i) {
        if (std::string(argv[i])=="-scale") scale = atof(argv[++i]);
		else if (std::string(argv[i])=="-width") w = atoi(argv[++i]);
		else if (std::string(argv[i])=="-height") h = atoi(argv[++i]);
    }
	return tracer::Pinhole(Eigen::Vector3f( 0.0f, 0.0f, -3.0f*scale), 
                        Eigen::Vector3f( 0.0f, 0.0f, 2.0f*scale), Eigen::Vector3f( 0, scale*float(h)/float(w), 0), Eigen::Vector3f(scale, 0, 0));
}

PointLight light_from_commandline(int argc, char** argv) {
    float scale = 1.0f; float energy = 1.0f;
    for (int i = 0; i<(argc-1);++i) {
        if (std::string(argv[i])=="-scale") scale = atof(argv[++i]);
    }
	for (int i = 0; i<argc;++i) {
		if (std::string(argv[i])=="arealight") energy=0.0f;
    }
	return PointLight(Eigen::Vector3f(1.5f*scale,0.0f,0.0f),
					SphericalSpectrumCone(Eigen::Vector3f(-1,0,0),M_PI/7.0f,
						Spectrum::Constant(energy)));
}

float error(const CImgWrapper<float>& gt, const CImgWrapper<float>& image) {
	float sum_num = 0.0f; float sum_den = 0.0f;
	for (auto p : multidimensional_range(gt.resolution())) {
		auto pgt = gt(p); auto pimage = image(p);
		for (int c = 0; c<3; ++c) {
			sum_num += (pgt[c] - pimage[c])*(pgt[c] - pimage[c]);
			sum_den += pgt[c]*pgt[c];
		}
	}
	return std::sqrt(sum_num)/std::sqrt(sum_den);
}

template<typename RF>
void compare_renderers(int argc, char** argv, const std::string& name, const RF& render_function, const CImgWrapper<float>& ground_truth) {
    int w = 512; int h = 512;
	float cv_rate = 0.5;
	unsigned long max_spp = 128;
	float error_rate = 1.e-5f;
    std::size_t seed = std::random_device()();
	
	for (int i = 0; i<(argc-1); ++i) {
		if (std::string(argv[i])=="-width") w = atoi(argv[++i]);
		else if (std::string(argv[i])=="-height") h = atoi(argv[++i]);
		else if (std::string(argv[i])=="-max-spp") max_spp = atol(argv[++i]);
		else if (std::string(argv[i])=="-cv-rate") cv_rate = atof(argv[++i]);
        else if (std::string(argv[i])=="-seed") seed = atol(argv[++i]);
		else if (std::string(argv[i])=="-error-rate") error_rate = atof(argv[++i]);
	}
	
	std::array<float,6> range_min, range_max; range_min.fill(0); range_max.fill(1);
	Range<float,6> render_range(range_min, range_max);
	CImgWrapper<float> image(w,h);

	std::cout<<name<<std::endl;
	for (unsigned long spp = 2; spp<2*max_spp; spp*=2) {
		image.clear();
		std::cout<<"  "<<std::setw(5)<<spp<<"spp "<<std::flush;
	    integrator_bins_stepper(stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),spp).integrate(image,image.resolution(),render_function, render_range);
		std::cerr<<"mc "<<error(ground_truth,image)<<" "<<std::flush;
		std::stringstream filename;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-montecarlo.hdr";
		image.save(filename.str());
		image.clear(); filename.str("");
		unsigned long spp_cv = std::min(spp-1,std::max(1UL,(unsigned long)(cv_rate*spp)));
		integrator_bins_stepper(stepper_bins_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(error_rate), (spp_cv*w*h)/(3*3*3*3*3*2), seed, seed+1),spp - spp_cv).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"cv "<<error(ground_truth,image)<<" "<<std::endl;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp.hdr";
		image.save(filename.str());
	}
}

int main(int argc, char** argv) {
	int w = 512; int h = 512;
	std::string output = "output";
	float scale = 1.0f;
	unsigned long gt_spp = 128000;
    std::size_t seed = std::random_device()();
	
	for (int i = 0; i<(argc-1); ++i) {
		if (std::string(argv[i])=="-width") w = atoi(argv[++i]);
		else if (std::string(argv[i])=="-height") h = atoi(argv[++i]);
		else if (std::string(argv[i])=="-ground-truth-spp") gt_spp = atol(argv[++i]);
		else if (std::string(argv[i])=="-scale") scale = atof(argv[++i]);
	    else if (std::string(argv[i])=="-seed") seed = atol(argv[++i]);
		else if (std::string(argv[i])=="-output") output = std::string(argv[++i]);
	}
	
	auto camera = camera_from_commandline(argc,argv);
	auto light = light_from_commandline(argc,argv);
	auto scene = scene_from_commandline(argc,argv);
	auto medium = medium_from_commandline(argc,argv);
	
	
	auto render_function = RenderMediumTwoBounces(scene,camera,medium,light,6.0f*scale);
	std::array<float,6> range_min, range_max; range_min.fill(0); range_max.fill(1);
	Range<float,6> render_range(range_min, range_max);
	
	std::string gt_filename = output+"-groundtruth-"+std::to_string(gt_spp)+"spp.hdr";
	CImgWrapper<float> ground_truth(w,h);
	
	if ((!ground_truth.load_hdr(gt_filename)) || 
			(std::array<std::size_t,2>{std::size_t(w),std::size_t(h)} != ground_truth.resolution())) {	
		integrate_bins_stepper_progression("Calculating ground truth ",stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),gt_spp, ground_truth,ground_truth.resolution(),render_function, render_range);
		ground_truth.save(gt_filename);
	} else std::cout<<"Ground truth loaded from "<<gt_filename<<std::endl;
	
	compare_renderers(argc,argv,output,render_function,ground_truth);
}




