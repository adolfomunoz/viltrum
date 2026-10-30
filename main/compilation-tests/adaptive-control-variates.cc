#include <iostream>
#include <iomanip>
#include "../../viltrum.h"
#include <cmath>
#include <string>

using namespace viltrum;
class DysonSeries {
    float decay;
public:
    DysonSeries(float decay = 0.75f) : decay(decay) {}
    template<typename Seq>
    float operator()(const Seq& seq) const {
        float sum = 1.0f; 
        float range = 1.0f;
        auto it = seq.begin();
        while ((*it) < decay) { //Multiply by decay, divide by decay as probability, gets simplified.
            ++it;
            float x = range*(*it);
            float prob = 1.0/range;
            range = x; //The new x reduces the range
            sum += x/prob;
            ++it;
        }
        return sum;
    }
} integrand_infinite;

template<typename Integrator, typename Function>
void test(const std::string& name, const Integrator& integrator, const Function& f, const std::size_t bins) {
    LoggerProgress logger(name);
    std::vector<float> sol(bins,0.0f);
    integrate(integrator,sol,f,range_primary_infinite<float>(),logger);
//    for (float v : sol) std::cout<<std::fixed<<std::setprecision(2)<<std::setw(4)<<v<<" ";
//    std::cout<<std::endl;
}


int main(int argc, char **argv) {

    std::size_t bins = 24;
    unsigned long iterations = 100;
    unsigned long mc_samples = 8;
    unsigned long spp_residual = 128;
    std::size_t seed = 1234;

    for (int i=1;i<argc;++i) {
        std::string arg = argv[i];
        if (arg == "--bins") { bins = std::stoul(argv[++i]); }
        else if (arg == "--iterations") { iterations = std::stoul(argv[++i]); }
        else if (arg == "--mc_samples") { mc_samples = std::stoul(argv[++i]); }
        else if (arg == "--spp_residual") { spp_residual = std::stoul(argv[++i]); }
        else if (arg == "--seed") { seed = std::stoul(argv[++i]); }
    }

    std::cout<<"Running adaptive control variates test with parameters: bins="<<bins<<", iterations="<<iterations<<", mc_samples="<<mc_samples<<", spp_residual="<<spp_residual<<", seed="<<seed<<"\n";
    test("Simple Monte Carlo",
        monte_carlo_per_bin(std::mt19937(seed),mc_samples*spp_residual),
        integrand_infinite,bins);


    test("Perbin stratification - control variates",
        integrator_adaptive_fubini_variance_reduction_parallel<4>(
            nested(simpson,trapezoidal),
            error_heuristic_size(error_metric_relative(),1.e-5),
            iterations,mc_samples,
            rr_uniform_region(), 
            cv_fixed_weight(1.0),
            region_sampling_uniform(),
            std::mt19937(seed),
            spp_residual),
        integrand_infinite,bins);
    
    test("Perregion stratification - control variates",
        integrator_adaptive_fubini_variance_reduction_parallel_regions<4>(
            nested(simpson,trapezoidal),
            error_heuristic_size(error_metric_relative(),1.e-5),
            iterations,mc_samples,
            region_stratification_uniform(),
            region_sampling_uniform(),
            std::mt19937(seed),
            spp_residual,1.0,1024),
        integrand_infinite,bins);

    test("Perregion stratification trapezoidal - control variates",
        integrator_adaptive_fubini_variance_reduction_parallel_regions_trapezoidal<4>(
            nested(simpson,trapezoidal),
            error_heuristic_size(error_metric_relative(),1.e-5),
            iterations,mc_samples,
            region_stratification_uniform(),
            region_sampling_uniform(),
            std::mt19937(seed),
            spp_residual,1.0,1024),
        integrand_infinite,bins);

	return 0;
}



