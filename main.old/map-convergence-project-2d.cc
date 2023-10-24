#include <memory>
#include "../plot/map-error-time.h"
#include "convergence-datafile.h"
#include "../quadrature/integrate-bins.h"
#include "../quadrature/integrate-bins-stepper.h"
#include "../quadrature/monte-carlo.h"
#include "../quadrature/control-variates.h"
#include "../plot/convergence.h"
#include "../functions/functions2d.h"
#include "../quadrature/integrate-adaptive-control-variates.h"
#include "../quadrature/integrate-optimized-adaptive-stratified-control-variates.h"

#include <iostream>
#include <cmath>
#include <chrono>

//Returns error-time
template<typename F, typename Error>
std::tuple<float,float> test(const F& f, const Error& error, 
	std::size_t mc_spp, std::size_t qdt_its, std::size_t bins, double error_size_weight) {

	std::vector<double> sol(bins,0.0); 	
	auto start = std::chrono::steady_clock::now();
	if (qdt_its == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_per_bin(stepper_monte_carlo_uniform()),
			std::max(mc_spp,std::size_t(1))
		),sol,f,range(0.0,0.0,1.0,1.0));
	} else if (mc_spp == 0) {
		integrate_bins(integrator_bins_stepper(
			stepper_bins_adaptive(nested(simpson,trapezoidal),
					error_single_dimension_size(error_size_weight)),
			qdt_its),sol,f,range(0.0,0.0,1.0,1.0));
	} else {
		integrate_bins(	integrator_optimized_perregion_adaptive_stratified_control_variates(
			nested(simpson, trapezoidal),
			error_single_dimension_size(error_size_weight),
			qdt_its, mc_spp),sol,f,range(0.0,0.0,1.0,1.0));	
	}
	auto end = std::chrono::steady_clock::now();
	return std::tuple<double,double>(
		error(sol),std::chrono::duration<double>(end-start).count());
}

int main(int argc, char **argv) {	
	std::string output = "output.svg";
	std::vector<int>  mc_spp{ 0, 4, 16, 64, 256 };
	std::vector<int> qdt_its{ 0, 2,  4,  8,  16 };
	auto [func,gt] = function2d(argc,argv);
//	double rmin = 0; double rmax = 1;
	double min_plottable_error = 1.e-10;
    std::size_t montecarlo_error_average = 100;
	std::size_t bins = 9;
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
		else if (std::string(argv[i])=="-bins") { bins = atoi(argv[++i]); }         
	}

	auto [ground_truth, persistent_data] = load_datafile(output.substr(0,output.length()-4)+".txt",bins);
	if (ground_truth.empty()) {
		ground_truth.resize(bins);
		std::cout<<"No data available nor ground truth available, generating"<<std::endl;
		double dx = 1.0/double(bins);
		for (std::size_t i=0; i<bins; ++i) ground_truth[i]=double(bins)*gt(dx*i,0,dx*(i+1),1);
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


	auto f = [func] (const std::array<double,2>& x) {
		return func(x[0],x[1]);
	};
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
					auto [e,t] = test(f,error,mc_spp[m],qdt_its[q],
									bins,error_size_weight);
					local_error+=e; local_time+=t;
				}
				time[m][q] = std::log(local_time/montecarlo_error_average);
				err[m][q] = std::log(local_error/montecarlo_error_average);
				persistent_data[mc_spp[m]][qdt_its[q]] = std::tuple<float,float>(time[m][q],err[m][q]);
			}
		}
		std::cout<<std::endl;
	}
		
		
	std::cout<<"Saving ground truth and data"<<std::endl;
	save_datafile(output.substr(0,output.length()-4)+".txt",
		ground_truth,persistent_data);
		
	std::ofstream file(output);
	file<<map_error_time(mc_spp,qdt_its,err,time,bins,"magma");
	return 0;
}



