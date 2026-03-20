#pragma once
#include <random>
#include "norm.h"
#include <vector>
#include "../range.h"

namespace viltrum {

class region_stratification_uniform {
public:
    template<typename Regions, typename RNG>
    std::vector<unsigned long> samples_per_region(
        unsigned long samples, const Regions& regions, RNG& rng) const {
        
        std::size_t nregions = regions.size();
        std::vector<unsigned long> spr(nregions,samples/nregions);
        unsigned long remainder = samples % nregions;
        if (remainder>0) {
            unsigned long random_start = std::uniform_int_distribution<unsigned long>(0, nregions - 1)(rng);
            for (unsigned long i=0;i<remainder;++i) 
                spr[(random_start+i)%nregions]++;
        }
        return spr;
    }
};




}   

