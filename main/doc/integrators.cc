#include <iostream>
#include <iomanip>
#include "../../viltrum.h"
#include <cmath>

class SquareHypotenuse {
public:
    float operator()(float x, float y) const {
        return x*x + y*y;
    }
} integrand;

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

int main(int argc, char **argv) {
    auto range = viltrum::range_all<2>(0.0f,1.0f);
    auto range_infinite = viltrum::range_primary_infinite<float>();
    std::vector<float> sol_bins(5,0.0f);

    //Monte-Carlo
    std::cout<<viltrum::integrate(viltrum::monte_carlo(std::ranlux48(),100), integrand, range)<<" ";
    std::cout<<viltrum::integrate(viltrum::monte_carlo(100,0), integrand_infinite,range_infinite)<<" ";
    viltrum::integrate(viltrum::monte_carlo(100), sol_bins, integrand, range);
    for (std::size_t i=0;i<sol_bins.size();++i) { std::cout<<" Bin "<<i<<": "<<sol_bins[i]<<" | ";}
    std::cout<<"\n";

    //Per-bins
    sol_bins.assign(sol_bins.size(),0.0f);
    viltrum::integrate(viltrum::integrator_per_bin(viltrum::monte_carlo(20,0)), sol_bins, integrand, range);
    for (std::size_t i=0;i<sol_bins.size();++i) { std::cout<<" Bin "<<i<<": "<<sol_bins[i]<<" | ";}
    sol_bins.assign(sol_bins.size(),0.0f);
    viltrum::integrate(viltrum::integrator_per_bin_parallel(viltrum::monte_carlo(20,0)), sol_bins, integrand_infinite, range_infinite);
    for (std::size_t i=0;i<sol_bins.size();++i) { std::cout<<" Bin "<<i<<": "<<sol_bins[i]<<" | ";}
    std::cout<<"\n";


    //Newton-Cotes Quadrature rules
    std::cout<<viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::trapezoidal), integrand, range)<<" ";
    std::cout<<viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::simpson), integrand, range);
    std::cout<<viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::boole), integrand, range)<<" ";
    sol_bins.assign(sol_bins.size(),0.0f);
    viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::steps<2>(viltrum::boole)), sol_bins, integrand, range);
    for (std::size_t i=0;i<sol_bins.size();++i) { std::cout<<" Bin "<<i<<": "<<sol_bins[i]<<" | ";}
    std::cout<<"\n";


    //Adaptive Nested Newton-Cotes with tolerance
    std::cout<<viltrum::integrate(viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),1.e-3), integrand, range)<<" "; 
    sol_bins.assign(sol_bins.size(),0.0f);
    viltrum::integrate(viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::boole,viltrum::simpson),viltrum::error_heuristic_size(viltrum::error_metric_relative(),1.e-5),1.e-3), sol_bins, integrand, range);
    for (std::size_t i=0;i<sol_bins.size();++i) { std::cout<<" Bin "<<i<<": "<<sol_bins[i]<<" | ";}
    std::cout<<"\n";

    //Adaptive Nested Newton-Cotes with number of steps
    std::cout<<viltrum::integrate(viltrum::integrator_adaptive_iterations(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),32), integrand, range)<<" "; 
    sol_bins.assign(sol_bins.size(),0.0f);
    viltrum::integrate(viltrum::integrator_adaptive_iterations_parallel(viltrum::nested(viltrum::boole,viltrum::simpson),viltrum::error_heuristic_size(viltrum::error_metric_relative(),1.e-5),32), sol_bins, integrand, range);
    for (std::size_t i=0;i<sol_bins.size();++i) { std::cout<<" Bin "<<i<<": "<<sol_bins[i]<<" | ";}
    std::cout<<"\n";

    //Fubini with adaptive newton cotes and Monte-Carlo
    std::cout<<viltrum::integrate(
        viltrum::integrator_fubini<1>(
            viltrum::integrator_adaptive_iterations_parallel(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),32),
            viltrum::monte_carlo(32)
        ), integrand, range)<<" ";

    //Fubini with adaptive newton cotes and Monte-Carlo, infinite dimensionality
    std::cout<<viltrum::integrate(
        viltrum::integrator_fubini<1>(
            viltrum::integrator_adaptive_iterations_parallel(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),32),
            viltrum::monte_carlo(32)
        ), integrand_infinite, range_infinite)<<"\n";

    
    //Crespo et al. 2021 integrator (with small improvements)
    std::cout<<viltrum::integrate(
        viltrum::integrator_crespo2021(16,64), integrand, range)<<" ";
    //Crespo et al. 2021 integrator for innite dimensionality (4 dimensions on the polynomial) 
    std::cout<<viltrum::integrate(
        viltrum::integrator_crespo2021_infinite<4>(16,4,64), integrand_infinite, range_infinite)<<"\n";

	return 0;
}



