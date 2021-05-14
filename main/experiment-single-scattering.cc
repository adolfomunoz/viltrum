#include "../functions/functions-render3d.h"
#include "../quadrature/range.h"
#include "../quadrature/integrate-adaptive-control-variates.h"
#include "../quadrature/integrate-bins-stepper.h"
#include "../quadrature/integrate-bins-adaptive.h"
#include "../quadrature/integrate-bins-adaptive-precalculate.h"
#include "../quadrature/integrate-optimized-adaptive-stratified-control-variates.h"
#include "../quadrature/munoz2014.h"
#include "../utils/cimg-wrapper.h"
#include "../plot/integration2d.h"
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <sstream>

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
    float scale = 1.0f; 
    for (int i = 0; i<(argc-1);++i) {
        if (std::string(argv[i])=="-scale") scale = atof(argv[++i]);
    }
	return PointLight(Eigen::Vector3f(1.5f*scale,0.0f,0.0f),
					SphericalSpectrumCone(Eigen::Vector3f(-1,0,0),M_PI/7.0f,
						Spectrum::Constant(1.0f)));
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

void add_labels(svg_cpp_plot::Graph2D& graph, int res) {
	graph.add(svg_cpp_plot::_2d::text({res/2,res+14},"u")).font_size(30).text_anchor(svg_cpp_plot::text_anchor_middle).dominant_baseline(svg_cpp_plot::dominant_baseline_hanging);
	graph.add(svg_cpp_plot::_2d::text({-14,res/2},"y")).font_size(30).text_anchor(svg_cpp_plot::text_anchor_end).dominant_baseline(svg_cpp_plot::dominant_baseline_middle);	
}

template<typename F>
svg_cpp_plot::SVG plot_cv_approximation(const F& f, float error_rate, unsigned long iterations, int res = 400) {
	svg_cpp_plot::SVG svg;
	auto& graph = svg.add(svg_cpp_plot::Graph2D({res,res},svg_cpp_plot::BoundingBox(0,0,1,1)));
	graph.area().add(plot_adaptive_approximation_2d(nested(simpson,trapezoidal),error_single_dimension_size(error_rate),f,range(0.0f,0.0f,1.0f,1.0f), iterations, res));
	graph.area().add(plot_adaptive_boundaries_2d(nested(simpson,trapezoidal),error_single_dimension_size(error_rate),f,range(0.0f,0.0f,1.0f,1.0f), iterations)).stroke_width(1).stroke(svg_cpp_plot::green);
	add_labels(graph,res);
	svg.viewBox(svg_cpp_plot::BoundingBox(-34,0,res,res+34));
	return svg;
}	

template<typename RF>
void compare_renderers(int argc, char** argv, const std::string& name, const RF& render_function, const CImgWrapper<float>& ground_truth) {
    int w = 512; int h = 512;
	float cv_rate = 0.5;
	unsigned long max_spp = 128;
	unsigned long spp_pixel = 4;
	float error_rate = 1.e-5f;
    std::size_t seed = std::random_device()();
	int slice_at_column = -1;
	
	for (int i = 0; i<(argc-1); ++i) {
		if (std::string(argv[i])=="-width") w = atoi(argv[++i]);
		else if (std::string(argv[i])=="-height") h = atoi(argv[++i]);
		else if (std::string(argv[i])=="-slice-at-column") slice_at_column = atoi(argv[++i]);
		else if (std::string(argv[i])=="-max-spp") max_spp = atol(argv[++i]);
		else if (std::string(argv[i])=="-cv-rate") cv_rate = atof(argv[++i]);
		else if (std::string(argv[i])=="-spp-pixel") spp_pixel = atol(argv[++i]);
        else if (std::string(argv[i])=="-seed") seed = atol(argv[++i]);
		else if (std::string(argv[i])=="-error-rate") error_rate = atof(argv[++i]);
	}
	
	if (slice_at_column < 0) slice_at_column = h/2;
	
	std::array<float,3> range_min, range_max; range_min.fill(0); range_max.fill(1);
	Range<float,3> render_range(range_min, range_max);
	CImgWrapper<float> image(w,h);
	
	
	float slice_max = 0.0f;

	for (int i = 0; i<h; ++i) for (int j = 0; j<h; ++j) {
		auto pixel = render_function(std::array<float,3>{float(i)/float(h), float(slice_at_column)/float(h),float(j)/float(h)});
		for (int c = 0; c<3; ++c) if (pixel[c]>slice_max) slice_max=pixel[c];
	}

	std::cerr<<"MAX = "<<slice_max<<std::endl;
	auto slice = [&] (const std::array<float,2>& x) -> Eigen::Array3f {
		return render_function(std::array<float,3>{x[0], float(slice_at_column)/float(h),x[1]})/slice_max;
	};

	std::cout<<name<<std::endl;
	for (unsigned long spp = 2; spp<2*max_spp; spp*=2) {
		
		image.clear();
		std::cout<<"  "<<std::setw(5)<<spp<<"spp "<<std::flush;
	    integrator_bins_stepper(stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),spp).integrate(image,image.resolution(),render_function, render_range);
		std::cerr<<"mc "<<error(ground_truth,image)<<" "<<std::flush;
		std::stringstream filename;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-montecarlo.hdr";
		image.save(filename.str());

		unsigned long spp_cv = std::min(spp-1,std::max(1UL,(unsigned long)(cv_rate*spp)));

		image.clear(); filename.str("");
		integrator_bins_stepper(stepper_bins_adaptive(nested(simpson,trapezoidal),error_single_dimension_size(error_rate)),(spp_cv*w*h)/(3*3*2)).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"ap "<<std::flush;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp-controlvariate";
		image.save(filename.str()+".hdr");
		std::ofstream(filename.str()+".svg") << plot_cv_approximation(slice,error_rate,(spp_cv*h)/(3*2),h);
			
		image.clear(); filename.str("");
		integrator_bins_stepper(stepper_bins_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(error_rate), (spp_cv*w*h)/(3*3*2), seed, seed+1),spp - spp_cv).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"cv1 "<<error(ground_truth,image)<<" "<<std::flush;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp-alpha1.hdr";
		image.save(filename.str());
					
		image.clear(); filename.str("");
		integrator_optimized_perpixel_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(error_rate), (spp_cv*w*h)/(3*3*2), spp - spp_cv, seed).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"cvOpp "<<error(ground_truth,image)<<" "<<std::flush;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp-alphaoptperpixel.hdr";
		image.save(filename.str());
		
		image.clear(); filename.str("");
		integrator_alpha1_perregion_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(error_rate), (spp_cv*w*h)/(3*3*2), spp - spp_cv, seed).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"cv1pr "<<error(ground_truth,image)<<" "<<std::flush;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp-alpha1perregion.hdr";
		image.save(filename.str());
							
		
		image.clear(); filename.str("");
		integrator_optimized_perregion_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(error_rate), (spp_cv*w*h)/(3*3*2), spp - spp_cv, seed).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"cvOpr "<<error(ground_truth,image)<<" "<<std::flush;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp-alphaoptperregion.hdr";
		image.save(filename.str());
							
		image.clear(); filename.str("");
		integrator_optimized_adaptive_stratified_control_variates(nested(simpson,trapezoidal),error_single_dimension_size(error_rate), (spp_cv*w*h)/(3*3*2), spp - spp_cv, seed).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"cvO "<<error(ground_truth,image)<<" "<<std::flush;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_cv<<"cvspp-alphaopt.hdr";
		image.save(filename.str());

		
		image.clear(); filename.str("");
        /*
		integrator_bins_stepper(stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),spp_pixel).integrate(image,image.resolution(),
			[&] (const std::array<float,2>& x) {
				return integrate(integrator_adaptive_iterations(nested(simpson,trapezoidal),error_single_dimension_size(error_rate),spp/spp_pixel),
				[&] (const std::array<float,1>& t) {
					return render_function(std::array<float,3>{x[0],x[1],t[0]});
				},range(0.0f,1.0f));
			}, range(0.0f,0.0f,1.0f,1.0f));
            */
        integrator_bins_munoz_2014(spp,spp_pixel, error_rate, seed).integrate(image,image.resolution(),render_function, render_range);
		std::cout<<"mu "<<error(ground_truth,image)<<" "<<std::endl;
		filename<<name<<"-"<<std::setfill('0')<<std::setw(5)<<spp<<"spp-"<<spp_pixel<<"perpixel-munoz2014.hdr";
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
	
	
	auto render_function_dumb = RenderMediumSingleScattering(scene,camera,medium,light,6.0f*scale);
	auto render_function_distance = RenderMediumSingleScatteringDistance(scene,camera,medium,light,6.0f*scale);
	auto render_function_equiangular = RenderMediumSingleScatteringEquiangular(scene,camera,medium,light,6.0f*scale);
	
	std::array<float,3> range_min, range_max; range_min.fill(0); range_max.fill(1);
	Range<float,3> render_range(range_min, range_max);
	
	std::string gt_filename = output+"-groundtruth-"+std::to_string(gt_spp)+"spp.hdr";
	CImgWrapper<float> ground_truth(w,h);
	
	if ((!ground_truth.load_hdr(gt_filename)) || 
			(std::array<std::size_t,2>{std::size_t(w),std::size_t(h)} != ground_truth.resolution())) {	
		integrate_bins_stepper_progression("Calculating ground truth ",stepper_bins_per_bin(stepper_monte_carlo_uniform(seed)),gt_spp, ground_truth,ground_truth.resolution(),render_function_dumb, render_range);
		ground_truth.save(gt_filename);
	} else std::cout<<"Ground truth loaded from "<<gt_filename<<std::endl;
	
	compare_renderers(argc,argv,output+"-equiangular",render_function_equiangular,ground_truth);
	compare_renderers(argc,argv,output+"-distance",render_function_distance,ground_truth);
	compare_renderers(argc,argv,output+"-dumb",render_function_dumb,ground_truth);

}




