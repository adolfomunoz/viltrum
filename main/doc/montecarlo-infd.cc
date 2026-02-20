#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

int main() {
    float decaying_factor = 0.75f;
    // This integrand is an infinite sum of decaying products.
    auto integrand = [decaying_factor] (const auto& seq) -> float {
        auto x = seq.begin(); //The sequence is an infinite sequence of numbers within the integration range, which is generated on the fly by the integrator. We can iterate over it until we want.
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
    };


    float sol = viltrum::integrate(
        viltrum::monte_carlo(8192),  //Numerical integration techinque: Monte-Carlo with 8192 samples
        integrand,
        viltrum::range_infinite<float>(0.0f,1.0f) //Range of infinite dimensionality, all dimensions between 0 and 1
    );

    // The integral of this function is a geometric series with ratio 0.5, so the result should be 2.0f.
    std::cout<<"Integral: "<<std::setprecision(6)<<sol<<" should be close to "<<(decaying_factor/(1.0f-decaying_factor))<<std::endl;
}