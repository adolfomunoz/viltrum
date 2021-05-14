#include "../multiarray/multiarray.h"
#include "../multiarray/transform.h"
#include "../multiarray/io.h"
#include <iostream>

float avg(const std::array<float,5>& p,std::size_t i) {
	if (i==0) return (p[0]+p[1])/2.0f;
    else if (i==4) return (p[3]+p[4])/2.0f;
    else return (p[i-1]+p[i]+p[i+1])/3.0f;
}

int main(int argc, char **argv) {	
	multiarray<float,5,2> ma(0);
	ma[{1,1}]=5; ma[{0,1}]=2; ma[{2,2}]=-10;
	std::cout<<ma<<std::endl;
	std::cout<<ma.transform(avg,0)<<std::endl;
	std::cout<<ma.transform(avg,1)<<std::endl;
	std::cout<<ma.transform(avg,0).transform(avg,1)<<std::endl;
	std::cout<<ma.transform_all(avg)<<std::endl;
	auto ma0 = clone(ma.transform(avg,0));
	std::cout<<ma0<<std::endl;
	auto ma01 = clone(ma.transform(avg,0).transform(avg,1));
	std::cout<<ma01<<std::endl;
	auto maall = clone(ma.transform_all(avg));
	std::cout<<maall<<std::endl;
	return 0;
}



