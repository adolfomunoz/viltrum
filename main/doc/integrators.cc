#include <iostream>
#include <iomanip>
#include "../viltrum.h"
#include <cmath>

float function(const std::array<float,2>& x) {
    return std::cos(4*M_PI*x[0]/(x[1]+1.0f));
}

int main(int argc, char **argv) {
    auto range = viltrum::range_all<2>(0.0f,1.0f);

    //Monte-Carlo
    std::cout<<viltrum::integrator_monte_carlo_uniform(std::ranlux48(),100).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_monte_carlo_uniform(100,0).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_monte_carlo_uniform(100).integrate(function,range)<<"\n";
    
    //Newton-Cotes Quadrature rules
    std::cout<<viltrum::integrator_quadrature(viltrum::trapezoidal).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_quadrature(viltrum::simpson).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_quadrature(viltrum::boole).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_quadrature(viltrum::steps<2>(viltrum::boole)).integrate(function,range)<<"\n";
    
    //Adaptive Nested Newton-Cotes with tolerance
    std::cout<<viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),1.e-3).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::boole,viltrum::simpson),1.e-3).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::steps<2>(viltrum::boole),viltrum::boole)).integrate(function,range)<<"\n";
 
    //Adaptive Nested Newton-Cotes with number of steps
    std::cout<<viltrum::integrator_adaptive_iterations(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),10).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_iterations(viltrum::nested(viltrum::boole,viltrum::simpson),10).integrate(function,range)<<"\n";

    //Adaptive Nested Newton-Cotes control variates
    std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),4,std::ranlux48(),80).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,std::ranlux48(),80).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),4,80,0).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,std::ranlux48(),80,0).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),4,80).integrate(function,range)<<" ";
    std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,std::ranlux48(),80).integrate(function,range)<<"\n";

	return 0;
}



