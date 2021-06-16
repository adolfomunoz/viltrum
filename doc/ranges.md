# `viltrum` - Integration range definition

Defining integration ranges in `viltrum` for numerical integration methods is relatively simple, with several options making it more flexible. However, their definition affects the whole integration because they not only define the range but the data type (float or double) that will be used to explored the integrand as well as the dimensionality of the problem. An integration range is defined as follows:their integration with any numerical method is simple. 

An integration range is represented by the class `viltrum::Range` that can be constructed as follows:

```cpp
viltrum::Range<F,N> integration_range(const std::array<F,N>& a, const std::array<F,N>& b)
``` 
where:
- `N` is the number of dimensions explored by the integration method, and should be the same number of dimensions than the [function that defines the integrand](integrands.md). It should be the same number than the dimensions of the integration range that is defined when integrating.
- `F` is a floating point value that explores the function (commonly `float` or `double`), and again must match the [integrand's parameter](integrands.md).
- `a` and `b` mark, respectively, the lower and upper bound of the integration range per dimension. 

For simplicity reasons, several helper functions have been defined to construct the integration range. The following examples will be using the following integrand:

```cpp
class Gaussian {
public:
    template<typename real, std::size_t NDIM>
    real operator()(const std::array<real,NDIM>& x) const {
        real sum(0);
        for (real xi : x) sum += xi*xi;
        return std::exp(-0.5*sum)/std::sqrt(2.0*M_PI);
    }
};
```

Note that the integrand above is defined for any number of dimensions `NDIM` as well as for any floating point data type `real` per dimension, so the integration range will define how this function is explored.

The simplest helper function `viltrum::range` requires two `std::arrays` as parameters, from which the  dimensionality and the real type for `viltrum::Range` are deduced. The following example uses this formulation for integration ranges of 1, 2 and 3 dimensions:
```cpp
int main(int argc, char** argv) {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::array<double,1> min1d{0.0}; std::array<double,1> max1d{1.0};
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range(min1d,max1d))<<std::endl;
    std::array<double,2> min2d{0.0,0.0}; std::array<double,2> max2d{1.0,1.0};
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range(min2d,max2d))<<std::endl;
    std::array<double,3> min3d{0.0,0.0,0.0}; std::array<double,3> max3d{1.0,1.0,1.0};
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range(min3d,max3d))<<std::endl;
}
```

Defining first the `std::array`s and then the integration range is quite boilerplate-y, so the `viltrum::range` function can be called with the per-dimesion limits directly (first the real numbers for the lower limits and then the real numbers for the upper limits) as illustrated in the following example:
```cpp
int main(int argc, char** argv) {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range(0.0, 1.0))<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range(0.0,0.0, 1.0,1.0))<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range(0.0,0.0,0.0, 1.0,1.0,1.0))<<std::endl;
}
```

Obviously, the number of floating point parameters must be even.

If the integration range is the same for all the dimensions, then it is possible to use the helper function `viltrum::range_all<NDIM>(a,b)`, where `NDIM` is the number of dimensions and `a` and `b` are floating point numbers that mark the lower and upper limits of the integration range for all the dimensions. This is illustrated in the following: 
```cpp
int main(int argc, char** argv) {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range_all<1>(0.0,1.0))<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range_all<2>(0.0,1.0))<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range_all<3>(0.0,1.0))<<std::endl;
}
```

Last helper function is `viltrum::range_primary<NDIM,real>()` which generates a `NDIM`-dimensional range of `real` floating point type between between 0 and 1 for all dimensions. If omitted, the `real` parameter defaults to `float`. It can be used as follows:
```cpp
int main(int argc, char** argv) {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),viltrum::range_primary<1>())<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),viltrum::range_primary<2>())<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),viltrum::range_primary<3>())<<std::endl;
}
```

The last option to generate ranges is to use the operator `|`, which concatenates two ranges, returning a higher dimensionality integration range, used as follows:
```cpp
int main(int argc, char** argv) {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    
    std::cout<<"1D = "<<method.integrate(Gaussian(),
            viltrum::range(0.0,1.0))<<std::endl;
    std::cout<<"2D = "<<method.integrate(Gaussian(),
            viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0))<<std::endl;
    std::cout<<"3D = "<<method.integrate(Gaussian(),
            viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0))<<std::endl;
}
```

The code illustrated in this page can be tested and compiled in a [source code example](../main/doc/ranges.cc)
 


