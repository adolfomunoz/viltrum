#include "../../viltrum.h"
#include <iostream>


class Series {
    float decaying_factor;
public:
    Series(float df = 0.5) : decaying_factor(df) {}

    //The parameter is a sequence of samples from the integration range, which is generated on the fly by the integrator. We can iterate over it until we want. In this case, we will iterate until we find a sample that is above the decaying factor, which is the Russian Roulette termination condition for this series.
    template<typename Sequence>
    float operator()(const Sequence& seq) const {
        auto x = seq.begin(); //This is the iterator over the sequence of samples.
        float sum = 0.0f;
        float term = 1.0f;  
        // For each iteration we take two samples:
        // - the first is Russian Roulette for the component with the decaying factor.
        // - the second is the value of the component with the decaying factor, which is multiplied by the term and added to the sum.       
        while ((*x) < decaying_factor) { 
            ++x; 
            term *= 2.0f*(*x); 
            ++x; 
            sum += term; 
        }
        return sum;
    }
};


int main(int argc, char** argv) {
    unsigned long samples = 16;
    for (int i = 0; i<argc-1; ++i) {
        if (std::string(argv[i])=="-samples") samples = atol(argv[++i]);
    }
    
    auto method = viltrum::monte_carlo(samples);
    std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_primary_infinite<float>())<<std::endl;
    
    std::vector<float> rangemin{0.0f,0.0f}, rangemax{1.0f,1.0f};
    std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_infinite(rangemin,rangemax))<<std::endl;
    
    std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_infinite<float>({0.0f},{1.0f}))<<std::endl;
    std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_infinite<float>({0.0f,0.0f},{1.0f,1.0f}))<<std::endl;

    std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range(0.0f,1.0f) | viltrum::range_primary_infinite<float>())<<std::endl;
}