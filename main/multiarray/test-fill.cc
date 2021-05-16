#include <iostream>
#include "../multiarray/multiarray.h"
#include "../multiarray/io.h"
#include "../multiarray/fold.h"
#include "../multiarray/fill.h"
#include <cmath>
#include <sstream>
#include <algorithm>

using namespace viltrum;

int main(int argc, char **argv) {
	auto f1 = [] (const std::array<double,1>& x) { return std::cos(M_PI*x[0]/2.0); };
	auto f2 = [] (const std::array<double,2>& x) { return x[0]*x[1]; };
	auto f3 = [] (const std::array<double,3>& x) { return x[0]*x[1]*x[2]; };
	multiarray<double,9,1> ma1;
	multiarray<double,5,2> ma2;
	multiarray<double,3,3> ma3;
	ma1.fill(f1); 
	ma2.fill(f2);
	ma3.fill(f3);
	
	std::cout<<ma1<<std::endl<<std::endl;
	std::cout<<ma2<<std::endl;
	std::cout<<ma3<<std::endl;
	
	ma3.fill(1);
	std::cout<<ma3<<std::endl;
	


	return 0;
}



