#include "../../src/newton-cotes/region.h"
#include "../../src/newton-cotes/region-trapezoidal.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace viltrum;

constexpr std::size_t NPOINTS = 1024;
volatile double sink = 0; //Avoids the compiler discarding or hoisting the benchmarked calls

//Average nanoseconds per call of f(i)
template<typename F>
double time_ns(std::size_t iterations, const F& f) {
    double acc = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i<iterations; ++i) acc += f(i);
    auto t1 = std::chrono::steady_clock::now();
    sink = sink + acc;
    return std::chrono::duration<double,std::nano>(t1-t0).count()/double(iterations);
}

template<typename R>
void benchmark(const std::string& name, const R& r, std::size_t iterations,
        const std::vector<std::array<float,R::dimensions>>& pos,
        const std::vector<std::array<float,R::dimensions>>& a,
        const std::vector<std::array<float,R::dimensions>>& b) {
    double t_app = time_ns(iterations, [&] (std::size_t i) { return double(r.approximation_at(pos[i%NPOINTS])); });
    double t_sub = time_ns(iterations, [&] (std::size_t i) { return double(r.integral_subrange(a[i%NPOINTS],b[i%NPOINTS])); });
    std::cout<<"  "<<std::left<<std::setw(24)<<name<<std::right
             <<std::setw(18)<<std::fixed<<std::setprecision(2)<<t_app
             <<std::setw(18)<<t_sub<<std::endl;
}

template<std::size_t D>
void run(std::size_t iterations) {
    std::array<float,D> xmin, xmax;
    xmin.fill(0); xmax.fill(1);
    auto f = [] (const std::array<float,D>& x) { 
        float s = 0; 
        for (float v : x) s += std::sin(3*v) + v*v; 
        return s; 
    };

    std::mt19937 rng(1);
    std::uniform_real_distribution<float> dist(0,1);
    std::vector<std::array<float,D>> pos(NPOINTS), a(NPOINTS), b(NPOINTS);
    for (std::size_t k = 0; k<NPOINTS; ++k)
        for (std::size_t d = 0; d<D; ++d) {
            pos[k][d] = dist(rng);
            a[k][d] = dist(rng);
            b[k][d] = std::min(1.0f,a[k][d] + 0.5f*dist(rng));
        }

    std::cout<<"Dimensions = "<<D<<std::endl;
    std::cout<<"  "<<std::left<<std::setw(24)<<"Region"<<std::right
             <<std::setw(18)<<"approximation_at"<<std::setw(18)<<"integral_subrange"<<" (ns/call)"<<std::endl;
    benchmark("Region trapezoidal", region(f,trapezoidal,xmin,xmax), iterations, pos, a, b);
    benchmark("Region simpson", region(f,simpson,xmin,xmax), iterations, pos, a, b);
    benchmark("Region boole", region(f,boole,xmin,xmax), iterations, pos, a, b);
    benchmark("RegionTrapezoidal", region_trapezoidal(f,xmin,xmax), iterations, pos, a, b);
}

int main(int argc, char** argv) {
    std::size_t iterations = 200000;
    for (int i = 0; i<(argc-1); ++i) 
        if (std::string(argv[i])=="-iterations") iterations = std::atol(argv[++i]);

    run<1>(iterations);
    run<2>(iterations);
    run<3>(iterations);
    run<4>(iterations);
    return 0;
}
