# `viltrum` - Integration range definition

## Finite dimensionality ranges
Defining finite integration ranges in `viltrum` for numerical integration methods is relatively simple, with several options making it more flexible. However, their definition affects the whole integration because they not only define the range but the data type (float or double) that will be used to explored the integrand as well as the dimensionality of the problem. An integration range is defined as follows:their integration with any numerical method is simple. 

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
auto method = viltrum::monte_carlo(1024);

std::array<double,1> min1d{0.0}; std::array<double,1> max1d{1.0};
std::cout<<"1D = "<<viltrum::integrate(method,Gaussian(),viltrum::range(min1d,max1d))<<std::endl;
std::array<double,2> min2d{0.0,0.0}; std::array<double,2> max2d{1.0,1.0};
std::cout<<"2D = "<<viltrum::integrate(method,Gaussian(),viltrum::range(min2d,max2d))<<std::endl;
std::array<double,3> min3d{0.0,0.0,0.0}; std::array<double,3> max3d{1.0,1.0,1.0};
std::cout<<"3D = "<<viltrum::integrate(method,Gaussian(),viltrum::range(min3d,max3d))<<std::endl;
```

Defining first the `std::array`s and then the integration range is quite boilerplate-y, so the `viltrum::range` function can be called with the per-dimesion limits directly (first the real numbers for the lower limits and then the real numbers for the upper limits) as illustrated in the following example:
```cpp
auto method = viltrum::monte_carlo(1024);

std::cout<<"1D = "<<viltrum::integrate(method,Gaussian(),viltrum::range(0.0, 1.0))<<std::endl;
std::cout<<"2D = "<<viltrum::integrate(method,Gaussian(),viltrum::range(0.0,0.0, 1.0,1.0))<<std::endl;
std::cout<<"3D = "<<viltrum::integrate(method,Gaussian(),viltrum::range(0.0,0.0,0.0, 1.0,1.0,1.0))<<std::endl;
```

Obviously, the number of floating point parameters must be even.

If the integration range is the same for all the dimensions, then it is possible to use the helper function `viltrum::range_all<NDIM>(a,b)`, where `NDIM` is the number of dimensions and `a` and `b` are floating point numbers that mark the lower and upper limits of the integration range for all the dimensions. This is illustrated in the following: 
```cpp
auto method = viltrum::monte_carlo(1024);

std::cout<<"1D = "<<viltrum::integrate(method,Gaussian(),viltrum::range_all<1>(0.0,1.0))<<std::endl;
std::cout<<"2D = "<<viltrum::integrate(method,Gaussian(),viltrum::range_all<2>(0.0,1.0))<<std::endl;
std::cout<<"3D = "<<viltrum::integrate(method,Gaussian(),viltrum::range_all<3>(0.0,1.0))<<std::endl;
```

Last helper function is `viltrum::range_primary<NDIM,real>()` which generates a `NDIM`-dimensional range of `real` floating point type between between 0 and 1 for all dimensions. If omitted, the `real` parameter defaults to `float`. It can be used as follows:

```cpp
auto method = viltrum::monte_carlo(1024);

std::cout<<"1D = "<<viltrum::integrate(method,Gaussian(),viltrum::range_primary<1>())<<std::endl;
std::cout<<"2D = "<<viltrum::integrate(method,Gaussian(),viltrum::range_primary<2>())<<std::endl;
std::cout<<"3D = "<<viltrum::integrate(method,Gaussian(),viltrum::range_primary<3>())<<std::endl;
```

The last option to generate ranges is to use the operator `|`, which concatenates two ranges, returning a higher dimensionality integration range, used as follows:
```cpp
auto method = viltrum::integrator_monte_carlo_uniform(1024);

std::cout<<"1D = "<<viltrum::integrate(method,Gaussian(),
        viltrum::range(0.0,1.0))<<std::endl;
std::cout<<"2D = "<<viltrum::integrate(method,Gaussian(),
        viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0))<<std::endl;
std::cout<<"3D = "<<viltrum::integrate(method,Gaussian(),
        viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0) | viltrum::range(0.0,1.0))<<std::endl;
```

The code illustrated here can be tested and compiled in a [source code example](../main/doc/ranges.cc)
 
## Infinite dimensionality ranges

Defining an infinite integration ranges in `viltrum` for numerical integration methods Their definition affects the integration process depending on the numerical data type (float or double) that will be used to explore the integration domain. Note that the integrand must be defined in a very particular way in order to support an infinite domain. Not all numerical integration methods support infinite integration ranges. An infinite integration range is defined as follows:

An integration range is represented by the class `viltrum::RangeInifnite` that can be constructed as follows:

```cpp
viltrum::RangeInfinite<F> integration_range(const std::vector<F>& a, const std::vector<F>& b)
``` 
where:
- `F` is a floating point value that explores the function (commonly `float` or `double`).
- `a` and `b` mark, respectively, the lower and upper bound of the integration range per dimension for the first set of dimensions. The rest of the dimensions, not present in the vectors, are always from 0 to 1.

The following examples will be using the following integrand, which fulfills the requirement of having an infinite sequence as parameter, by being a template:

```cpp
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
```

The most common infinite-dimensionality domain is the primary domain, with the limits for all dimensions being between 0 and 1:

```cpp
auto method = viltrum::monte_carlo(samples);
std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_primary_infinite<float>())<<std::endl;
```

However, the standard way of defining the range is through two vectors. As the vectors are of limited size, the limits for all dimensions outside of the vector are, once again, 0 and 1.

```cpp
auto method = viltrum::monte_carlo(samples);
std::vector<float> rangemin{0.0f,0.0f}, rangemax{1.0f,1.0f};
std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_infinite(rangemin,rangemax))<<std::endl;
```

These vectors can also be defined inline, for which the range needs to know the datatype in advance, and cannot deduce it:

```cpp
auto method = viltrum::monte_carlo(samples);
std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_infinite<float>({0.0f},{1.0f}))<<std::endl;
std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range_infinite<float>({0.0f,0.0f},{1.0f,1.0f}))<<std::endl;
```

The last option to generate ranges is to use the operator `|`, which concatenates two ranges, and the rightmost one can be an infinite range, returning another infinite range with the limits combined:

```cpp
auto method = viltrum::monte_carlo(samples);
std::cout<<"Integral = "<<viltrum::integrate(method,Series(),viltrum::range(0.0f,1.0f) | viltrum::range_primary_infinite<float>())<<std::endl;
```

The code illustrated here can be tested and compiled in a [source code example](../main/doc/ranges-infinite.cc)
