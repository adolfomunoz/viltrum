# `viltrum` - Integrand definition

Defining integrands in `viltrum` for their integration with any numerical method is simple. You need to define a function that has the form:

```cpp
V f(const std::array<F,N>& x)
``` 
where:
- `N` is the number of dimensions explored by the function. It should be the same number than the dimensions of the (integration range)[ranges.md] that is defined when integrating.
- `F` is a floating point value that explores the function (commonly `float` or `double`) and must match the data type of the (integration range)[ranges.md].
- `V` is the resulting value of the function. `V` needs to have some numeric operations: addition `+`, multiplication by a scalar `*`, division by a scalar `/`... In general, floating point numbers and standard numerical algebra arrays (`Eigen::Array`) comply with this requirements.

There are several ways in which such an integrand can be defined in C++, and each of them will be illustrated with an working example. All of the examples below approximate the volume of a sphere using Monte Carlo integration.

Integrands can be defined: 
- Through a standard C++ function:

```cpp
#include <viltrum/viltrum.h>
float sphere(const std::array<float,3>& x) {
    return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<=1.0f?1.0f:0.0f;
}
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate(sphere,range)<<std::endl;
}
```

- Through a class that represents a function (defining its `operator()`):

```cpp
#include <viltrum/viltrum.h>
class Sphere {
public:
    float operator()(const std::array<float,3>& x) const {
        return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<=1.0f?1.0f:0.0f;        
    }
};
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate(Sphere(),range)<<std::endl;
}
```

- Through a lambda expression:

```cpp
#include <viltrum/viltrum.h>
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate([] (const std::array<float,3>& x) { return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1.0f?1.0f:0.0f; },range)<<std::endl;
}
```

## Wrapping normal functions

It is true that defining functions having `std::array`s as parameters is rather unusual and unintuitive. More often you will define a function as having multiple independent parameters. In order to fit those as `viltrum` integrands, the library provides a wrapper that automatically does that, which is unsurprisingly named `viltrum::function_wrapper(...)`. You can wrap any function in the invocation of the integration method.

Lets see how to apply the wrapper with similar examples than above:

- With a C++ function:
```cpp
#include <viltrum/viltrum.h>
float sphere_parameters(float x, float y, float z) {
    return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f;
}
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate(viltrum::function_wrapper(sphere_parameters),range)<<std::endl;
}
```

- With a C++ class representing a function with an `operator()`:
```cpp
#include <viltrum/viltrum.h>
class SphereParameters {
public:
    float operator()(float x, float y, float z) const {
    return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f;
}
};
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate(viltrum::function_wrapper(SphereParameters()),range)<<std::endl;
}
```

- With a lambda expression:
```cpp
#include <viltrum/viltrum.h>
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate(viltrum::function_wrapper([] (float x, float y, float z) { return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f; }),
                        range)<<std::endl;
}
```

The code illustrated in this page can be tested and compiled in a [source code example](../main/doc/integrands.cc)


