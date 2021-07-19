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
  
/*
	std::cout<<std::endl<<"Adaptive tolerance"<<std::endl;
	test("Tr.Simp. e-2  ",integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-2),f,rmin,rmax);
	test("Tr.Simp. e-4  ",integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-4),f,rmin,rmax);
	test("Tr.Simp. e-6  ",integrator_adaptive_tolerance(nested(simpson,trapezoidal),1.e-6),f,rmin,rmax);
	test("Simp.Boole e-2",integrator_adaptive_tolerance(nested(boole,simpson),1.e-2),f,rmin,rmax);
	test("Simp.Boole e-4",integrator_adaptive_tolerance(nested(boole,simpson),1.e-4),f,rmin,rmax);
	test("Simp.Boole e-6",integrator_adaptive_tolerance(nested(boole,simpson),1.e-6),f,rmin,rmax);

	std::cout<<std::endl<<"Adaptive iterations"<<std::endl;
	test("Tr.Simp.    1 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),1),f,rmin,rmax);
	test("Tr.Simp.   10 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),10),f,rmin,rmax);
	test("Tr.Simp.  100 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),100),f,rmin,rmax);
	test("Tr.Simp. 1000 ",integrator_adaptive_iterations(nested(simpson,trapezoidal),1000),f,rmin,rmax);
	test("Simp.Boole   1",integrator_adaptive_iterations(nested(boole,simpson),1),f,rmin,rmax);
	test("Simp.Boole  10",integrator_adaptive_iterations(nested(boole,simpson),10),f,rmin,rmax);
	test("Simp.Boole 100",integrator_adaptive_iterations(nested(boole,simpson),100),f,rmin,rmax);
	test("Simp.Boole1000",integrator_adaptive_iterations(nested(boole,simpson),1000),f,rmin,rmax);


	std::cout<<std::endl<<"Monte Carlo"<<std::endl;
	test("Uniform     1 ",integrator_monte_carlo_uniform(1),f,rmin,rmax);
	test("Uniform    10 ",integrator_monte_carlo_uniform(10),f,rmin,rmax);
	test("Uniform   100 ",integrator_monte_carlo_uniform(100),f,rmin,rmax);
	test("Uniform  1000 ",integrator_monte_carlo_uniform(1000),f,rmin,rmax);
*/

	return 0;
}



