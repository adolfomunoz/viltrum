#include <iostream>
#include "../multiarray/multiarray.h"
#include "../multiarray/io.h"
#include "../multiarray/split.h"
#include <cmath>
#include <sstream>
#include <algorithm>

template<std::size_t size>
void test1(std::size_t parts) {
	unsigned int evals = 0;
	auto f = [&evals] (const std::array<double,1>& x) { ++evals; return x[0]; };
	multiarray<double,size,1> m; m.fill(f);
	std::cout<<m<<" -> \t";
	for (auto p : m.split(f,0,parts)) std::cout<<"| "<<p<<" |\t";
	std::cout<<"Evaluations = "<<evals<<std::endl;
}

template<std::size_t size>
void test2(std::size_t parts) {
	unsigned int evals = 0;
	auto f = [&evals] (const std::array<double,2>& x) { ++evals; return x[0]*x[1]; };
	multiarray<double,size,2> m; m.fill(f);
	std::cout<<std::endl<<m<<std::endl;
	for (auto p : m.split(f,0,parts)) std::cout<<p<<"----------"<<std::endl;
	std::cout<<"Evaluations = "<<evals<<std::endl<<std::endl;
}

int main(int argc, char **argv) {
	test1<3>(2);
	test1<5>(2);
	test1<3>(3);
	test2<3>(2);
	test2<5>(2);
	test2<3>(3);
	return 0;
}



