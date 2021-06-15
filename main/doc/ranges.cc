#include "../../viltrum.h"


class Gaussian {
public:
    template<typename real, std::size_t NDIM>
    real operator()(const std::array<real,NDIM>& x) const {
        real sum(0);
        for (real xi : x) sum += xi*xi;
        return std::exp(-0.5*sum)/std::sqrt(2.0*M_PI);
    }
};

int main(int argc, char** argv) {
    unsigned long samples = 16;
    for (int i = 0; i<argc-1; ++i) {
        if (std::string(argv[i])=="-samples") samples = atol(argv[++i]);
    }
    
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::array<double,1> min1d{0.0}; std::array<double,1> max1d{1.0};
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range(min1d,max1d))<<std::endl;
    std::array<double,2> min2d{0.0,0.0}; std::array<double,2> max2d{1.0,1.0};
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range(min2d,max2d))<<std::endl;
    std::array<double,3> min3d{0.0,0.0,0.0}; std::array<double,3> max3d{1.0,1.0,1.0};
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range(min3d,max3d))<<std::endl;
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range(0.0,1.0))<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range(0.0,0.0,1.0,1.0))<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range(0.0,0.0,0.0,1.0,1.0,1.0))<<std::endl;
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range(0.0,1.0))<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),
            viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0))<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),
            viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0))<<std::endl;
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range_all<1>(0.0,1.0))<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range_all<2>(0.0,1.0))<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range_all<3>(0.0,1.0))<<std::endl;
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range_primary<1>())<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range_primary<2>())<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range_primary<3>())<<std::endl;
}
