#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

int main() {
    float sol[10]; //We will compute the integral in 10 bins, so we need an array of 10 floats to store the results
    //Access to the binning structure (sol) through the position (pos, the index) in the range. 
    auto sol_access = [&sol] (const std::array<std::size_t,1>& pos) -> float& { return sol[pos[0]]; }; 
    //  The function to integrate is f(x,y) = x^2 + y^2. The parameter is an std::array with the two dimenions but it could also be a bidimensional function
    auto integrand = [] (const std::array<float,2>& x) -> float { return x[0]*x[0] + x[1]*x[1]; };
    //Integration range: from (0,0) to (1,1). Could also be viltrum::range(0.0f,1.0f,0.0f,1.0f)
    auto range = viltrum::range(std::array<float,2>{0.0f,0.0f}, std::array<float,2>{1.0f,1.0f}); 
    viltrum::integrate(
        viltrum::monte_carlo(8192),  //Numerical integration techinque: Monte-Carlo with 8192 samples
        sol_access,
        std::array<std::size_t,1>{10}, //Number of bins in each dimension (8 bins in this case) integrand. Dimensionality should be the same than the parameter of the access above
        integrand,
        range
    );

    for (std::size_t i=0;i<10;++i) { std::cout<<"Bin "<<i<<": "<<sol[i]<<"\n";}
}