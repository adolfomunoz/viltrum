#pragma once
#include <random>
#include "norm.h"
#include <vector>
#include "../range.h"

namespace viltrum {

class rr_uniform_region {
public:
    class RR {
    private:
        //This is the MonteCarlo residual, RR among regions
        std::uniform_int_distribution<std::size_t> rr;
        template<typename Rs>
        RR(const Rs& regions) : rr(0,regions.size()-1) {}
    public:
        template<typename RNG>
        std::tuple<std::size_t,double> choose(RNG& rng) { 
            return std::tuple(rr(rng),rr.max()+1); 
        }
        friend class rr_uniform_region;
    };
    template<typename Regions>
    RR russian_roulette(const Regions& regions) const {
        return RR(regions);
    }
};

template<typename Norm = NormDefault>
class rr_integral_region {
    Norm norm;
public:
    rr_integral_region(const Norm& n = Norm()) : norm(n) {}
    class RR {
    private:
        //This is the MonteCarlo residual, RR among regions
        std::discrete_distribution<std::size_t> rr;
        std::vector<double> weights;
        Norm norm;

        template<typename Rs>
        RR(const Rs& regions, const Norm& n) : weights(regions.size()), norm(n) {
            std::size_t i = 0; 
            for (const auto& [r,region_bin_range] : regions) {
                weights[i++] = norm(r->integral_subrange(region_bin_range));
            }
            //This is to avoid zeroes and negative numbers: put a minimum relative weight
            double sum = 0.0;
            for (double w : weights) sum += w;
            if (sum<=0.0) for(double& w : weights) w = 1.0;
            else for (double& w : weights) w = std::max(w,0.01*sum/double(weights.size()));
            rr = std::discrete_distribution<std::size_t>(weights.begin(),weights.end());
        }
        
    public:
        template<typename RNG>
        std::tuple<std::size_t,double> choose(RNG& rng) { 
            std::size_t choice = rr(rng);
            return std::tuple(choice,1.0/rr.probabilities()[choice]); 
        }
        friend class rr_integral_region;
    };
    template<typename Regions>
    RR russian_roulette(const Regions& regions) const {
        return RR(regions,norm);
    }
};

template<typename Norm = NormDefault>
class rr_error_region {
    Norm norm;
public:
    rr_error_region(const Norm& n = Norm()) : norm(n) {}
    class RR {
    private:
        //This is the MonteCarlo residual, RR among regions
        std::discrete_distribution<std::size_t> rr;
        std::vector<double> weights;
        Norm norm;

        template<typename Rs>
        RR(const Rs& regions, const Norm& n) : weights(regions.size()), norm(n) {
            std::size_t i = 0; 
            for (const auto& [r,region_bin_range] : regions) {
                weights[i++] = norm(r->error())*region_bin_range.volume()/r->range().volume();
            }
            //This is to avoid zeroes and negative numbers: put a minimum relative weight
            double sum = 0.0;
            for (double w : weights) sum += w;
            if (sum<=0.0) for(double& w : weights) w = 1.0;
            else for (double& w : weights) w = std::max(w,0.01*sum/double(weights.size()));
            rr = std::discrete_distribution<std::size_t>(weights.begin(),weights.end());
        }
        
    public:
        template<typename RNG>
        std::tuple<std::size_t,double> choose(RNG& rng) { 
            std::size_t choice = rr(rng);
            return std::tuple(choice,1.0/rr.probabilities()[choice]); 
        }
        friend class rr_error_region;
    };
    template<typename Regions>
    RR russian_roulette(const Regions& regions) const {
        return RR(regions,norm);
    }
};

template<typename Norm = NormDefault>
class rr_pdf_region {
    Norm norm;
public:
    rr_pdf_region(const Norm& n = Norm()) : norm(n) {}
    class RR {
    private:
        //This is the MonteCarlo residual, RR among regions
        std::discrete_distribution<std::size_t> rr;
        std::vector<double> weights;
        Norm norm;

        template<typename Rs>
        RR(const Rs& regions, const Norm& n) : weights(regions.size()), norm(n) {
            std::size_t i = 0; 
            for (const auto& [r,region_bin_range] : regions) {
                weights[i++] = r->pdf_integral_subrange(region_bin_range,norm);
            }
            //This is to avoid zeroes and negative numbers: put a minimum relative weight
            double sum = 0.0;
            for (double w : weights) sum += w;
            if (sum<=0.0) for(double& w : weights) w = 1.0;
            else for (double& w : weights) w = std::max(w,0.01*sum/double(weights.size()));
            rr = std::discrete_distribution<std::size_t>(weights.begin(),weights.end());
        }
        
    public:
        template<typename RNG>
        std::tuple<std::size_t,double> choose(RNG& rng) { 
            std::size_t choice = rr(rng);
            return std::tuple(choice,1.0/rr.probabilities()[choice]); 
        }
        friend class rr_pdf_region;
    };
    template<typename Regions>
    RR russian_roulette(const Regions& regions) const {
        return RR(regions,norm);
    }
};




}   

