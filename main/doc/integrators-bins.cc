#include <iostream>
#include <iomanip>
#include "../viltrum.h"
#include <cmath>

float slope(const std::array<float,2>& x) { return (x[1]<x[0])?1.0f:0.0f; }

int main(int argc, char **argv) {
    auto range = viltrum::range(std::array<float,2>{0,0},std::array<float,2>{1,1});
    std::vector<float> output(16);
    
    //Monte-Carlo
    viltrum::integrate_bins(
        viltrum::integrator_bins_monte_carlo_uniform(std::ranlux48(),100), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_bins_monte_carlo_uniform(100,0), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_bins_monte_carlo_uniform(100), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;
    std::cout<<std::endl;
    
    
    //Adaptive (only iterations in this case)
    viltrum::integrate_bins(
        viltrum::integrator_bins_adaptive(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),
        viltrum::error_relative_single_dimension(),100), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_bins_adaptive(viltrum::nested(viltrum::boole,viltrum::simpson),100), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    std::cout<<std::endl;
    

    //Per bin
    viltrum::integrate_bins(
        viltrum::integrator_bins_per_bin(viltrum::integrator_quadrature(viltrum::trapezoidal)), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;    
    viltrum::integrate_bins(
        viltrum::integrator_bins_per_bin(viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::boole,viltrum::simpson),1.e-3)), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;    
    viltrum::integrate_bins(
        viltrum::integrator_bins_per_bin(viltrum::integrator_monte_carlo_uniform(100)), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    std::cout<<std::endl;  

    //Optimized adaptive stratified control variates
    viltrum::integrate_bins(
        viltrum::integrator_optimized_perpixel_adaptive_stratified_control_variates(
            viltrum::nested(viltrum::simpson,viltrum::trapezoidal), // nested rule, order 2 polynomials
            viltrum::error_single_dimension_size(1.e-5), // error heuristic
            16, // number of adaptive iterations calculated from the spps
            100, // number of spps for the residual
            std::ranlux48()            // random number generator
        ), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_optimized_perpixel_adaptive_stratified_control_variates(
            viltrum::nested(viltrum::simpson,viltrum::trapezoidal), // nested rule, order 2 polynomials
            viltrum::error_single_dimension_size(1.e-5), // error heuristic
            16, // number of adaptive iterations calculated from the spps
            100, // number of spps for the residual
            0 // seed (can be omitted)
        ), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_alpha1_perpixel_adaptive_stratified_control_variates(
            viltrum::nested(viltrum::simpson,viltrum::trapezoidal), // nested rule, order 2 polynomials
            viltrum::error_single_dimension_size(1.e-5), // error heuristic
            16, // number of adaptive iterations calculated from the spps
            100 // seed (can be omitted)
        ), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_optimized_perregion_adaptive_stratified_control_variates(
            viltrum::nested(viltrum::simpson,viltrum::trapezoidal), // nested rule, order 2 polynomials
            viltrum::error_single_dimension_size(1.e-5), // error heuristic
            16, // number of adaptive iterations calculated from the spps
            100, // number of spps for the residual
            0 // seed (can be omitted)
        ), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    viltrum::integrate_bins(
        viltrum::integrator_alpha1_perregion_adaptive_stratified_control_variates(
            viltrum::nested(viltrum::simpson,viltrum::trapezoidal), // nested rule, order 2 polynomials
            viltrum::error_single_dimension_size(1.e-5), // error heuristic
            16, // number of adaptive iterations calculated from the spps
            100 // seed (can be omitted)
        ), 
        output, slope, range);
    for (float f : output) std::cout<<f<<" ";  
    std::cout<<std::endl;  
    std::cout<<std::endl;  

	return 0;
}



