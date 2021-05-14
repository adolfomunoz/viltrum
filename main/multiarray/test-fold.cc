#include "../multiarray/multiarray.h"
#include "../multiarray/fold.h"
#include "../multiarray/io.h"
#include "../quadrature/rules.h"
#include <iostream>

float sum(const std::array<float,3>& p) {
	return p[0]+p[1]+p[2];
}

int main(int argc, char **argv) {	
	multiarray<float,3,2> ma(0);
	ma[{1,1}]=5; ma[{0,1}]=2; ma[{2,2}]=-10;
	std::cout<<ma.fold(sum)<<std::endl;
	std::cout<<ma.fold(sum).fold(sum)<<std::endl;
	std::cout<<ma.fold_all(sum)<<std::endl;
	
	auto constant_function = [] (const std::array<double,3>& p)
	{	return 1; };
	
	multiarray<double,3,3> fevals;
	for (std::size_t i=0;i<3;++i) for (std::size_t j=0;j<3;++j) for (std::size_t k=0;k<3;++k) {
		fevals[{i,j,k}] = constant_function({i/2.0,j/2.0,k/2.0});
	}
	std::cout<<"Simpson = "<<fevals.fold_all(simpson)<<std::endl;
	std::cout<<"Simpson at 0.3,0.3,0.3 = "<<
		fevals.fold([&] (const auto& v) { return simpson.at(0.3,v);})
		      .fold([&] (const auto& v) { return simpson.at(0.3,v);})
		      .fold([&] (const auto& v) { return simpson.at(0.3,v);})
		      <<std::endl;
	
	return 0;
}



