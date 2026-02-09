#include "../../viltrum.h"
#include <iostream>
#include <iomanip>

int main() {
    float sol = viltrum::integrate(
        viltrum::monte_carlo(64),  //Numerical integration techinque: Monte-Carlo with 64 samples
        [] (float x) -> float { return std::sin(x);}, //Function to integrate: sin(x)
        viltrum::range(0.0f,3.14159265f) //Integration range: from 0 to pi
    );

    std::cout<<std::fixed<<std::setprecision(6)<<sol<<" should be close to 2\n";
}