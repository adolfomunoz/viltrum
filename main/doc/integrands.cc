#include "../../viltrum.h"

float sphere(const std::array<float,3>& x) {
    return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<=1.0f?1.0f:0.0f;
}

class Sphere {
public:
    float operator()(const std::array<float,3>& x) const {
        return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<=1.0f?1.0f:0.0f;        
    }
};

float sphere_parameters(float x, float y, float z) {
    return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f;
}

class SphereParameters {
public:
    float operator()(float x, float y, float z) const {
    return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f;
}
};

int main(int argc, char** argv) {
    unsigned long samples = 1024;
    for (int i = 0; i<argc-1; ++i) {
        if (std::string(argv[i])=="-samples") samples = atol(argv[++i]);
    }
    
    auto method = viltrum::monte_carlo(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    
    std::cout<<"Sphere volume = "<<(4.0f*M_PI/3.0f)<<std::endl;
    std::cout<<viltrum::integrate(method,sphere,range)<<std::endl;
    std::cout<<viltrum::integrate(method,Sphere(),range)<<std::endl;
    std::cout<<viltrum::integrate(method,[] (const std::array<float,3>& x) { return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1.0f?1.0f:0.0f; },range)<<std::endl;

    std::cout<<viltrum::integrate(method,sphere_parameters,range)<<std::endl;
    std::cout<<viltrum::integrate(method,SphereParameters(),range)<<std::endl;
    std::cout<<viltrum::integrate(method,[] (float x, float y, float z) { return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f; },range)<<std::endl;
}
