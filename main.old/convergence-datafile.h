#pragma once

#include <fstream>

//Ground truth and data. We mark the lack of a file as an empty ground truth
std::tuple< std::vector<double>,
   std::map<std::size_t,std::map<std::size_t,std::tuple<float,float>>> > load_datafile(std::string_view filename, std::size_t bins = 1) {
	std::map<std::size_t,std::map<std::size_t,std::tuple<float,float>>> persistent_data;
	
	std::ifstream fdatain; 
	fdatain.open(std::string(filename));
	if (!fdatain.good()) 
		return std::make_tuple(std::vector<double>(0),persistent_data);
	else {
		std::vector<double> ground_truth(bins);
		{ //Loading ground truth
			std::string line_str;
			std::getline(fdatain,line_str);
			std::stringstream line(line_str);
			for (auto& gt : ground_truth)
				line>>gt;
		}
		while (!fdatain.eof()) {
			std::string line_str;
			std::getline(fdatain,line_str);
			std::stringstream line(line_str);
			std::size_t mc, qdt; float log_time, log_error;
			line >> mc >> qdt >> log_time >> log_error;
			persistent_data[mc][qdt] = std::tuple<float,float>(log_time, log_error);
		}
		return std::make_tuple(ground_truth,persistent_data);
	}
}

void save_datafile(std::string_view filename,
	const std::vector<double>& ground_truth,
	const std::map<std::size_t,std::map<std::size_t,std::tuple<float,float>>>& persistent_data) {
		
	std::ofstream fdataout;
	fdataout.open(std::string(filename));
	if (fdataout.good()) {
		for (auto gtv : ground_truth)
			fdataout<<gtv<<" ";
		fdataout<<std::endl;
		for (auto [mc,inner] : persistent_data)
			for (auto [qdt,t] : inner) {
				fdataout<<mc<<" "<<qdt<<" "<<std::get<0>(t)<<" "<<std::get<1>(t)<<std::endl;
			}
	}
}