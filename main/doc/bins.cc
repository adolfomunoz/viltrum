#include <iostream>
#include <iomanip>
#include "../viltrum.h"
#include <cmath>


float slope(const std::array<float,2>& x) { return (x[1]<x[0])?1.0f:0.0f; }

int main(int argc, char **argv) {
    auto integrator_bins = viltrum::integrator_bins_monte_carlo_uniform(1024);
    auto range = viltrum::range(std::array<float,2>{0,0},std::array<float,2>{1,1});
    float output_array[16];
    
    auto output_array_access = [&output_array] (const std::array<std::size_t,1>& i) -> float& { return output_array[i[0]]; }; 
    integrator_bins.integrate(output_array_access,std::array<std::size_t,1>{16},slope,range);
    for (float f : output_array) std::cout<<f<<" ";
    std::cout<<std::endl;
    
    integrate_bins(integrator_bins,output_array_access,std::array<std::size_t,1>{16},slope,range);
    for (float f : output_array) std::cout<<f<<" ";
    std::cout<<std::endl;
    
    std::vector<float> output_vector(16);
    integrate_bins(integrator_bins, output_vector, slope, range);
    for (float f : output_vector) std::cout<<f<<" ";
    std::cout<<std::endl;
    
    std::vector<std::vector<float>> output_matrix(16,std::vector<float>(16));
    integrate_bins(integrator_bins, output_matrix, slope, range);
    std::cout<<std::endl;
    for (const std::vector<float>& row : output_matrix) {
        for (float f : row) std::cout<<f<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    
    
	return 0;
}



