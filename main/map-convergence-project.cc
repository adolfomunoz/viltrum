#include <memory>
#include "../plot/map-error-time.h"
#include "convergence-datafile.h"
#include "../quadrature/integrate-bins.h"
#include "../quadrature/integrate-bins-stepper.h"
#include "../quadrature/monte-carlo.h"
#include "../quadrature/control-variates.h"
#include "../plot/convergence.h"
#include "../functions/functions.h"
#include "../quadrature/vector-dimensions.h"
#include "../quadrature/integrate-adaptive-control-variates.h"
#include "../quadrature/integrate-optimized-adaptive-stratified-control-variates.h"

#include <iostream>
#include <cmath>
#include <chrono>

std::vector<double> generate_ground_truth(const Function1D& func, const std::array<std::size_t,4>& bins) {
    auto gt = std::get<1>(func);
    vector_dimensions<double,1> vd(std::array<std::size_t,1>{bins[0]},0.0);
    float dx = 1.0/float(bins[0]);
    for (std::size_t i = 0; i<bins[0]; ++i)
        vd(std::array{i}) = gt(i*dx,(i+1)*dx);
    return vd.raw_data();
}

std::vector<double> generate_ground_truth(const Function2D& func, const std::array<std::size_t,4>& bins) {
    auto gt = std::get<1>(func);
    vector_dimensions<double,2> vd(std::array<std::size_t,2>{bins[0],bins[1]},0.0);
    float dx = 1.0/float(bins[0]);
    float dy = 1.0/float(bins[1]);
    for (std::size_t i = 0; i<bins[0]; ++i)
        for (std::size_t j = 0; j<bins[1]; ++j)
            vd(std::array{i,j}) = gt(i*dx,j*dy,(i+1)*dx,(j+1)*dy);
    return vd.raw_data();
}

std::vector<double> generate_ground_truth(const Function3D& func, const std::array<std::size_t,4>& bins) {
    auto gt = std::get<1>(func);
    vector_dimensions<double,3> vd(std::array<std::size_t,3>{bins[0],bins[1],bins[2]},0.0);
    float dx = 1.0/float(bins[0]);
    float dy = 1.0/float(bins[1]);
    float dz = 1.0/float(bins[2]);
    for (std::size_t i = 0; i<bins[0]; ++i)
        for (std::size_t j = 0; j<bins[1]; ++j)
            for (std::size_t k = 0; k<bins[2]; ++k)
                vd(std::array{i,j,k}) = gt(i*dx,j*dy,k*dz,(i+1)*dx,(j+1)*dy,(k+1)*dz);
    return vd.raw_data();
}

std::vector<double> generate_ground_truth(const Function4D& func, const std::array<std::size_t,4>& bins) {
    auto gt = std::get<1>(func);
    vector_dimensions<double,4> vd(bins,0.0);
    float dx = 1.0/float(bins[0]);
    float dy = 1.0/float(bins[1]);
    float dz = 1.0/float(bins[2]);
	float dw = 1.0/float(bins[3]);
    for (std::size_t i = 0; i<bins[0]; ++i)
        for (std::size_t j = 0; j<bins[1]; ++j)
            for (std::size_t k = 0; k<bins[2]; ++k)
				for (std::size_t l = 0; l<bins[3]; ++l)
                vd(std::array{i,j,k,l}) = gt(i*dx,j*dy,k*dz,l*dw,(i+1)*dx,(j+1)*dy,(k+1)*dz,(l+1)*dw);
    return vd.raw_data();
}

//Returns error-time
template<typename Error>
std::tuple<float,float> test(const Function1D& func, const Error& error, 
	std::size_t mc_spp, std::size_t qdt_its, std::array<std::size_t,4> bins, double error_size_weight) {
    vector_dimensions<double,1> vd(std::array<std::size_t,1>{bins[0]},0.0);
    auto f = [func] (const std::array<double,1>& x) {
		return std::get<0>(func)(x[0]);
	};
	auto start = std::chrono::steady_clock::now();
	if (qdt_its == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_per_bin(stepper_monte_carlo_uniform()),
			std::max(mc_spp,std::size_t(1))
		),vd,vd.resolution(),f,range(0.0,1.0));
	} else if (mc_spp == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_adaptive(nested(simpson,trapezoidal),
					error_single_dimension_size(error_size_weight)),
			qdt_its),vd,vd.resolution(),f,range(0.0,1.0));
	} else {
		integrate_bins(	
		//integrator_optimized_perregion_adaptive_stratified_control_variates(
		integrator_optimized_perpixel_adaptive_stratified_control_variates(
			nested(simpson, trapezoidal),
			error_single_dimension_size(error_size_weight),
			qdt_its, mc_spp),vd,vd.resolution(),f,range(0.0,1.0));
	}
	auto end = std::chrono::steady_clock::now();
	return std::tuple<float,float>(
		error(vd.raw_data()),std::chrono::duration<float>(end-start).count());
}

//Returns error-time
template<typename Error>
std::tuple<float,float> test(const Function2D& func, const Error& error, 
	std::size_t mc_spp, std::size_t qdt_its, std::array<std::size_t,4> bins, double error_size_weight) {
    vector_dimensions<double,2> vd(std::array<std::size_t,2>{bins[0],bins[1]},0.0);
    auto f = [func] (const std::array<double,2>& x) {
		return std::get<0>(func)(x[0],x[1]);
	};
	auto start = std::chrono::steady_clock::now();
	if (qdt_its == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_per_bin(stepper_monte_carlo_uniform()),
			std::max(mc_spp,std::size_t(1))
		),vd,vd.resolution(),f,range(0.0,0.0,1.0,1.0));
	} else if (mc_spp == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_adaptive(nested(simpson,trapezoidal),
					error_single_dimension_size(error_size_weight)),
			qdt_its),vd,vd.resolution(),f,range(0.0,0.0,1.0,1.0));
	} else {
		integrate_bins(	
		//integrator_optimized_perregion_adaptive_stratified_control_variates(
		integrator_optimized_perpixel_adaptive_stratified_control_variates(
			nested(simpson, trapezoidal),
			error_single_dimension_size(error_size_weight),
			qdt_its, mc_spp),vd,vd.resolution(),f,range(0.0,0.0,1.0,1.0));
	}
	auto end = std::chrono::steady_clock::now();
	return std::tuple<float,float>(
		error(vd.raw_data()),std::chrono::duration<float>(end-start).count());
}

template<typename Error>
std::tuple<float,float> test(const Function3D& func, const Error& error, 
	std::size_t mc_spp, std::size_t qdt_its, std::array<std::size_t,4> bins, double error_size_weight) {
    vector_dimensions<double,3> vd(std::array<std::size_t,3>{bins[0],bins[1],bins[2]},0.0);
    auto f = [func] (const std::array<double,3>& x) {
		return std::get<0>(func)(x[0],x[1],x[2]);
	};
	auto start = std::chrono::steady_clock::now();
	if (qdt_its == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_per_bin(stepper_monte_carlo_uniform()),
			std::max(mc_spp,std::size_t(1))
		),vd,vd.resolution(),f,range(0.0,0.0,0.0,1.0,1.0,1.0));
	} else if (mc_spp == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_adaptive(nested(simpson,trapezoidal),
					error_single_dimension_size(error_size_weight)),
			qdt_its),vd,vd.resolution(),f,range(0.0,0.0,0.0,1.0,1.0,1.0));
	} else {
		integrate_bins(	
		//integrator_optimized_perregion_adaptive_stratified_control_variates(
		integrator_optimized_perpixel_adaptive_stratified_control_variates(
			nested(simpson, trapezoidal),
			error_single_dimension_size(error_size_weight),
			qdt_its, mc_spp),vd,vd.resolution(),f,range(0.0,0.0,0.0,1.0,1.0,1.0));
	}
	auto end = std::chrono::steady_clock::now();
	return std::tuple<float,float>(
		error(vd.raw_data()),std::chrono::duration<float>(end-start).count());
}

template<typename Error>
std::tuple<float,float> test(const Function4D& func, const Error& error, 
	std::size_t mc_spp, std::size_t qdt_its, std::array<std::size_t,4> bins, double error_size_weight) {
    vector_dimensions<double,4> vd(std::array<std::size_t,4>{bins[0],bins[1], bins[2],bins[3]},0.0);

    auto f = [func] (const std::array<double,4>& x) -> double {
		//return std::get<0>(func)(x[0],x[1],x[2],x[3]);
		double val = std::get<0>(func)(x[0],x[1],x[2],x[3]);
		//std::cout << "Val " << val << std::endl;
		return val;
	};
	auto start = std::chrono::steady_clock::now();
	if (qdt_its == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_per_bin(stepper_monte_carlo_uniform()),
			std::max(mc_spp,std::size_t(1))
		),vd,vd.resolution(),f,range(0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0));
	} else if (mc_spp == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_adaptive(nested(simpson,trapezoidal),
					error_single_dimension_size(error_size_weight)),
			qdt_its),vd,vd.resolution(),f,range(0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0));
	} else {
		integrate_bins(	
		//integrator_optimized_perregion_adaptive_stratified_control_variates(
		integrator_optimized_perpixel_adaptive_stratified_control_variates(
			nested(simpson, trapezoidal),
			error_single_dimension_size(error_size_weight),
			qdt_its, mc_spp),vd,vd.resolution(),f,range(0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0));
	}
	auto end = std::chrono::steady_clock::now();
	//std::cout << vd.raw_data()[0] << std::endl;
	return std::tuple<float,float>(
		error(vd.raw_data()),std::chrono::duration<float>(end-start).count());
}


int main(int argc, char **argv) {	
	std::string output = "output.svg";
	std::vector<int>  mc_spp{ 0, 4, 16, 64, 256 };
	std::vector<int> qdt_its{ 0, 2,  4,  8,  16 };
	auto func = function_from_commandline(argc,argv);
//	double rmin = 0; double rmax = 1;
	double min_plottable_error = 1.e-10;
    std::size_t montecarlo_error_average = 100;
	std::array<std::size_t,4> bins{9,1,1,1};
    double error_size_weight = 1.e-3;

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-montecarlo-spp") {
			std::cerr<<"MonteCarlo spp = ";
			mc_spp.clear(); ++i;
			for (;(i<argc)&&(argv[i][0]!='-');++i) {
				mc_spp.push_back(atol(argv[i]));
				std::cerr<<mc_spp.back()<<" ";
			}
			std::cerr<<std::endl; --i;
		}	
		else if (std::string(argv[i])=="-cv-iterations") {
			std::cerr<<"CV iterations = ";
			qdt_its.clear(); ++i;
			for (;(i<argc)&&(argv[i][0]!='-');++i) {
				qdt_its.push_back(atol(argv[i]));
				std::cerr<<qdt_its.back()<<" ";
			}
			std::cerr<<std::endl; --i;		
		}				
		else if (std::string(argv[i])=="-min-plottable-error") { min_plottable_error = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-montecarlo-error-average") { montecarlo_error_average = atof(argv[++i]); } 
		else if (std::string(argv[i])=="-error-size-weight") { error_size_weight=atof(argv[++i]); }
		else if (std::string(argv[i])=="-bins") { 
        	std::cerr<<"Bins = "; std::size_t b = 0;
            ++i;
			for (;(b<3)&&(i<argc)&&(argv[i][0]!='-');++i,++b) {
				bins[b]=atol(argv[i]);
				std::cerr<<bins[b]<<" ";
			}
			std::cerr<<std::endl; --i;	
        }         
	}

	auto [ground_truth, persistent_data] = load_datafile(output.substr(0,output.length()-4)+".txt",bins[0]*bins[1]*bins[2]*bins[3]);
	if (ground_truth.empty()) {
		std::cout<<"No data available nor ground truth available, generating"<<std::endl;
        ground_truth = std::visit([bins] (const auto& f) { return generate_ground_truth(f,bins); }, func);
	} else { std::cout<<"Loading data and ground truth"<<std::endl; }

    auto error = [ground_truth,min_plottable_error] (const std::vector<double>& a) {
        double err = 0;
        for (std::size_t i=0;i<a.size();++i)
            err += std::abs(ground_truth[i]-a[i]);
		return std::max(min_plottable_error,err/double(a.size()));
	};

    std::cerr<<"Sampling error average = "<<montecarlo_error_average<<std::endl;
	
	std::vector<std::vector<float>> time(mc_spp.size(),std::vector<float>(qdt_its.size()));
	std::vector<std::vector<float>> err(mc_spp.size(),std::vector<float>(qdt_its.size()));



	for (std::size_t m = 0; m<mc_spp.size(); ++m) {
		std::cout<<std::setw(6)<<mc_spp[m]<<" | "<<std::flush;
		for (std::size_t q = 0; q<qdt_its.size(); ++q) {
			std::cout<<std::setw(6)<<qdt_its[q]<<std::flush;
			//We just copy the value if it is already precalculated
			if ( (persistent_data.count(mc_spp[m]) > 0) &&
			   (persistent_data[mc_spp[m]].count(qdt_its[q]) > 0) ) { 
				std::tie(time[m][q],err[m][q]) = persistent_data[mc_spp[m]][qdt_its[q]];
				std::cout<<"^"<<std::flush;
			} else {
				std::cout<<"v"<<std::flush;
				float local_error(0), local_time(0);
				for (std::size_t a = 0; a<montecarlo_error_average; ++a) {
					auto [e,t] = 
                        std::visit([&] (const auto& f) { return test(f,error,mc_spp[m],qdt_its[q],bins,error_size_weight);},func);
					local_error+=e; local_time+=t;
				}
				time[m][q] = std::log(local_time/double(montecarlo_error_average));
				err[m][q] = std::log(local_error/double(montecarlo_error_average));
				persistent_data[mc_spp[m]][qdt_its[q]] = std::tuple<float,float>(time[m][q],err[m][q]);
			}
		}
		std::cout<<std::endl;
	}
		
		
	std::cout<<"Saving ground truth and data"<<std::endl;
	save_datafile(output.substr(0,output.length()-4)+".txt",
		ground_truth,persistent_data);
		
	std::ofstream file(output);
	file<<map_error_time(mc_spp,qdt_its,err,time,bins[0]*bins[1]*bins[2]*bins[3],"magma",function_dimensions(func));
	return 0;
}



