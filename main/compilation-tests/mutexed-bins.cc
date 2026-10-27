// Compares the execution time of MutexedBins and MutexedBinsAligned when several threads 
// accumulate random values into 2D bins.
// Usage: mutexed-bins [--threads N] [--locks N] [--adds N] [--resolution N] [--repetitions N]
// Run with --help for the description and defaults of each parameter.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <thread>
#include <vector>
#include "../../src/newton-cotes/mutexed-bins.h"
#include "../../src/newton-cotes/mutexed-bins-aligned.h"

class Bins2D {
    std::size_t width;
    std::vector<double> data;
public:
    Bins2D(std::size_t w, std::size_t h) : width(w), data(w*h,0.0) {}
    double& operator()(const std::array<std::size_t,2>& i) { return data[i[0]+width*i[1]]; }
    double sum() const { double s = 0; for (double d : data) s += d; return s; }
};

// Returns the elapsed time in seconds. Every add accumulates 1, so the sum of all bins must be equal to adds.
template<template<typename,std::size_t> class MutexedBinsType>
double run(std::size_t nthreads, std::size_t nlocks, std::size_t nadds, std::size_t resolution) {
    Bins2D bins(resolution,resolution);
    MutexedBinsType<Bins2D,2> mutexedbins(bins,{resolution,resolution},nlocks);
    auto start = std::chrono::steady_clock::now();
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < nthreads; ++t) {
        std::size_t count = nadds/nthreads + (t < nadds%nthreads ? 1 : 0);
        threads.emplace_back([&mutexedbins,t,count,resolution] {
            // Simple LCG so that random number generation does not dominate the measurement
            std::uint64_t s = 0x9E3779B97F4A7C15ULL*(t+1);
            for (std::size_t k = 0; k < count; ++k) {
                s = s*6364136223846793005ULL + 1442695040888963407ULL;
                mutexedbins.add({std::size_t((s>>20)%resolution),std::size_t((s>>40)%resolution)},1.0);
            }
        });
    }
    for (auto& th : threads) th.join();
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if (bins.sum() != double(nadds)) std::cerr<<"WARNING: wrong result, expected "<<nadds<<" and got "<<bins.sum()<<std::endl;
    return elapsed;
}

struct Options {
    std::size_t threads = std::max(1u,std::thread::hardware_concurrency());
    std::size_t locks = 16;
    std::size_t adds = 10000000;
    std::size_t resolution = 1024;
    std::size_t repetitions = 3;
};

void print_usage(const char* program) {
    Options d;
    std::cerr<<"Usage: "<<program<<" [--threads N] [--locks N] [--adds N] [--resolution N] [--repetitions N]\n"
             <<"  --threads      number of threads (default: "<<d.threads<<", the hardware concurrency)\n"
             <<"  --locks        number of locks (default: "<<d.locks<<")\n"
             <<"  --adds         total number of additions, split among the threads (default: "<<d.adds<<")\n"
             <<"  --resolution   the bins are resolution x resolution (default: "<<d.resolution<<")\n"
             <<"  --repetitions  number of times each version is run (default: "<<d.repetitions<<")\n"
             <<"Parameters can also be given as --name=value."<<std::endl;
}

// Returns false (after printing the reason) if the command line is not valid
bool parse_options(int argc, char** argv, Options& options) {
    const std::map<std::string,std::size_t Options::*> parameters = {
        {"threads",&Options::threads}, {"locks",&Options::locks}, {"adds",&Options::adds},
        {"resolution",&Options::resolution}, {"repetitions",&Options::repetitions} };
    for (int i = 1; i < argc; ++i) {
        std::string name = argv[i], value;
        if (name.rfind("--",0) != 0) { std::cerr<<"Unexpected argument '"<<name<<"'"<<std::endl; return false; }
        name = name.substr(2);
        auto eq = name.find('=');
        if (eq != std::string::npos) { value = name.substr(eq+1); name = name.substr(0,eq); }
        else if (i+1 < argc) value = argv[++i];
        auto it = parameters.find(name);
        if (it == parameters.end()) { std::cerr<<"Unknown parameter '--"<<name<<"'"<<std::endl; return false; }
        if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos) {
            std::cerr<<"Parameter '--"<<name<<"' requires a non-negative integer value"<<std::endl; return false;
        }
        options.*(it->second) = std::stoull(value);
    }
    if (options.threads == 0 || options.locks == 0 || options.resolution == 0 || options.repetitions == 0) {
        std::cerr<<"threads, locks, resolution and repetitions must be greater than zero"<<std::endl;
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) 
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") { print_usage(argv[0]); return 0; }
    Options options;
    if (!parse_options(argc,argv,options)) { print_usage(argv[0]); return 1; }
    const std::size_t nthreads = options.threads, nlocks = options.locks, nadds = options.adds, 
                      resolution = options.resolution, repetitions = options.repetitions;

    // Alternating the order of the versions avoids favouring one of them (warm-up, frequency scaling...)
    double plain = 0, aligned = 0;
    for (std::size_t r = 0; r < repetitions; ++r) {
        if (r%2 == 0) { 
            plain += run<viltrum::MutexedBins>(nthreads,nlocks,nadds,resolution);
            aligned += run<viltrum::MutexedBinsAligned>(nthreads,nlocks,nadds,resolution);
        } else {
            aligned += run<viltrum::MutexedBinsAligned>(nthreads,nlocks,nadds,resolution);
            plain += run<viltrum::MutexedBins>(nthreads,nlocks,nadds,resolution);
        }
    }
    plain /= repetitions; aligned /= repetitions;

    std::cout<<nthreads<<" threads, "<<nlocks<<" locks, "<<nadds<<" adds, "<<resolution<<"x"<<resolution<<" bins (mean of "<<repetitions<<" runs)"<<std::endl;
    std::cout<<"MutexedBins:        "<<plain<<" s"<<std::endl;
    std::cout<<"MutexedBinsAligned: "<<aligned<<" s"<<std::endl;
    std::cout<<"Speedup:            "<<plain/aligned<<"x"<<std::endl;
    return 0;
}
