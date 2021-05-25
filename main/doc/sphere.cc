#include "../../viltrum.h"


float sphere(const std::array<float,3>& x) {
    return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1.0f?1.0f:0.0f;
}

class Sphere {
public:
    float operator()(const std::array<float,3>& x) const {
        return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1.0f?1.0f:0.0f;        
    }
};



int main(int argc, char** argv) {
    unsigned long samples = 16;
    for (int i = 0; i<argc-1; ++i) {
        if (std::string(argv[i])=="-samples") samples = atol(argv[++i]);
    }
    
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::cout<<"Sphere volume = "<<(4.0f*M_PI/3.0f)<<" ";
    std::cout<<method.integrate(sphere,viltrum::range_all<3>(-1.0f,1.0f))<<" ";
    std::cout<<method.integrate(Sphere(),viltrum::range_all<3>(-1.0f,1.0f))<<" ";
    std::cout<<method.integrate([] (const std::array<float,3>& x) { return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1.0f?1.0f:0.0f; },
                    viltrum::range_all<3>(-1.0f,1.0f)))<<std::endl;
}