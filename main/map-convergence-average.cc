#include <memory>
#include "../plot/map-error-time.h"
#include "convergence-datafile.h"

#include <iostream>
#include <cmath>
#include <chrono>

int main(int argc, char **argv) {	
	std::string output = "output.svg";
	std::vector<int>  mc_spp{ 0, 4, 16, 64, 256 };
	std::vector<int> qdt_its{ 0, 2,  4,  8,  16 };
	std::vector<std::string> datafiles;
 	std::size_t bins = 9;
	bool geometric = true;

	for (int i = 0; i<argc; ++i) {
		if (std::string(argv[i])=="-geometric") geometric=true;
		else if (std::string(argv[i])=="-arithmetic") geometric=false;
	}

	for (int i = 0; i<argc-1; ++i) {
		if (std::string(argv[i])=="-output") { output = argv[++i]; } 
		else if (std::string(argv[i])=="-datafiles") {
			std::cerr<<"Data files = ";
			datafiles.clear(); ++i;
			for (;(i<argc)&&(argv[i][0]!='-');++i) {
				datafiles.push_back(std::string(argv[i]));
				std::cerr<<datafiles.back()<<" ";
			}
			std::cerr<<std::endl; --i;
		}	
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
		else if (std::string(argv[i])=="-bins") { bins = atoi(argv[++i]); }         
		
	}

	
 	std::vector<std::vector<float>> time(mc_spp.size(),std::vector<float>(qdt_its.size(),0));
	std::vector<std::vector<float>> err(mc_spp.size(),std::vector<float>(qdt_its.size(),0));
	std::vector<std::vector<int>> number(mc_spp.size(),std::vector<int>(qdt_its.size(),0));

	for (const auto& file : datafiles) {
		auto persistent_data = std::get<1>(load_datafile(file));
		for (std::size_t m = 0; m<mc_spp.size(); ++m) {
			for (std::size_t q = 0; q<qdt_its.size(); ++q) {
				if ( (persistent_data.count(mc_spp[m]) > 0) &&
			   (persistent_data[mc_spp[m]].count(qdt_its[q]) > 0) ) {
					auto [t,e] = persistent_data[mc_spp[m]][	qdt_its[q]];
					if (geometric) {
						time[m][q]+=t; 
						err[m][q]+=e; 
					} else {
						time[m][q]+=std::exp(t); 
						err[m][q]+=std::exp(e); 
					}
					number[m][q]++;
			   }
			}
		}
	}
	
	for (std::size_t m = 0; m<mc_spp.size(); ++m) {
		for (std::size_t q = 0; q<qdt_its.size(); ++q) {
			if (number[m][q] > 0) {
				if (geometric) {
					time[m][q]=time[m][q]/number[m][q];
					err[m][q]=err[m][q]/number[m][q];
				} else {
					time[m][q]=std::log(time[m][q]/number[m][q]);
					err[m][q]=std::log(err[m][q]/number[m][q]);
				}
			}
		}
	}
			
	std::ofstream file(output);
	file<<map_error_time(mc_spp,qdt_its,err,time,bins,"magma");
	return 0;
}



