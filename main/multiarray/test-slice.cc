#include "../multiarray/multiarray.h"
#include "../multiarray/io.h"
#include "../multiarray/slice.h"
#include "../multiarray/fill.h"
#include <iostream>

int main(int argc, char **argv) {	
	multiarray<float,3,2> m(0);
	std::cout<<m<<std::endl;
	m[{1,1}]=5;
	std::cout<<m<<std::endl;
	m.slice(0,0) = multiarray<float,3,1>(1);
	std::cout<<m<<std::endl;
	std::cout<<m.slice(1,1)<<std::endl<<std::endl;
	m.slice(2,1).fill(-1);
	std::cout<<m<<std::endl;

	return 0;
}



