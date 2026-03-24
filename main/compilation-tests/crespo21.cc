#include <iostream>
#include <iomanip>
#include "../../viltrum.h"
#include <cmath>

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
    auto range_infinite = viltrum::range_primary_infinite<float>();
    std::vector<float> sol_bins(10,0.0f);

    viltrum::integrate(viltrum::integrator_crespo2021_infinite<4>(1280,8,256), sol_bins, integrand_infinite, range_infinite);

    for (float x : sol_bins) std::cout<<x<<" ";
    std::cout<<std::endl;
	return 0;
}



