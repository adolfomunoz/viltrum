#include "../../src/monte-carlo/random-sequence-ref-dis.h"
#include "../../src/monte-carlo/random-sequence-ref.h"
#include "../../src/monte-carlo/random-sequence-rng.h"
#include "../../src/monte-carlo/random-sequence-seed.h"
#include "../../src/rng/pcg_random.hpp"
#include "../../src/rng/XoshiroCpp.hpp"
#include <iostream>
#include <chrono>
#include <string>
#include <cstdio>



using namespace viltrum;


template<template <typename, typename> class Seq>
struct sequence {
    template<typename RNG>
    static auto create(RNG& rng, unsigned seed, const RangeInfinite<float>& range) {
        return Seq<RNG,float>(rng, range);
    }
};

template<>
struct sequence<RandomSequenceSeed> {
    template<typename RNG>
    static auto create(RNG& rng, unsigned seed, const RangeInfinite<float>& range) {
        return RandomSequenceSeed<RNG,float>(seed, range);
    }
};

template<template <typename, typename> class Seq, typename RNG = std::mt19937>
void test_sequence(const char* name, const RangeInfinite<float>& range, unsigned seed = 0, unsigned long n = 10000000, unsigned long seqs = 10000) {
    RNG rng(seed);

    auto t0 = std::chrono::high_resolution_clock::now();
    unsigned long i = 0;
    float sum = 0.0f;
    {
        auto seq = sequence<Seq>::create(rng, seed, range);
        for (float x : seq) {
            sum += x;
            if (++i >= n) break;
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout << name<<"  - Single: "<<sum << "  (elapsed: " << elapsed << " s)";

    t0 = std::chrono::high_resolution_clock::now();
    {
        sum = 0.0f;
        for (unsigned long s = 0; s<seqs; ++s) {
            i = 0;
            auto seq = sequence<Seq>::create(rng, unsigned(rng()), range);
            for (float x : seq) {
                sum += x;
                if (++i >= (n / seqs)) break;
            }
        }
    }
    t1 = std::chrono::high_resolution_clock::now();

    elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cout <<"  - Multiple: "<<sum << "  (elapsed: " << elapsed << " s)" << std::endl;
}

int main(int argc, char** argv) {
    unsigned long nitems = 10000000;
    unsigned long nseqs  = 100000;

    for (int i = 0; i < argc-1; ++i) {
        if (std::string(argv[i])=="-nitems")     nitems = std::atol(argv[++i]);
        else if (std::string(argv[i])=="-nseqs") nseqs = std::atol(argv[++i]);
    }

    std::cout<<"NItems = "<<nitems<<" - NSeqs = "<<nseqs<<std::endl;

    RangeInfinite<float> range(std::vector{-2.0f},std::vector{2.0f});
    test_sequence<RandomSequenceSeed, std::mt19937>("Seed - mt19937", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRNG, std::mt19937>("RNG  - mt19937",  range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRef, std::mt19937>("Ref  - mt19937",  range,0, nitems, nseqs);
    test_sequence<RandomSequenceRefDis, std::mt19937>("Dis  - mt19937",  range,0, nitems, nseqs);
    
    test_sequence<RandomSequenceSeed, XoshiroCpp::Xoshiro128PlusPlus>("Seed - xoshi. ", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRNG, XoshiroCpp::Xoshiro128PlusPlus>("RNG  - xoshi. ", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRef, XoshiroCpp::Xoshiro128PlusPlus>("Ref  - xoshi. ", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRefDis, XoshiroCpp::Xoshiro128PlusPlus>("Dis  - xoshi. ", range, 0, nitems, nseqs);

    test_sequence<RandomSequenceSeed, pcg32>("Seed - pcg32   ", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRNG, pcg32>("RNG  - pcg32   ", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRef, pcg32>("Ref  - pcg32   ", range, 0, nitems, nseqs);
    test_sequence<RandomSequenceRefDis, pcg32>("Dis  - pcg32   ", range, 0, nitems, nseqs);
    std::cout<<std::endl;
}  