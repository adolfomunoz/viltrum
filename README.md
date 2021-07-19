# viltrum
VILTRUM: Varied Integration Layouts for arbiTRary integrals in a Unified Manner - A C++17 header-only library that provides a set of numerical integration routines.

## Installation

`viltrum` is a header-only library and has no external dependencies, so installation is rather simple (it does not require any compilation process). You just need to download the library in a specific folder:

```
folder> git clone https://github.com/adolfomunoz/viltrum.git
```

Then make sure `folder` is within the included directories (`-Ifolder` parameter in g++, `include_directories("folder")` in CMake) and 
include it from C++.

```cpp
#include <viltrum/viltrum.h>
```

That would be enough to use all the features of the library. There are other alternatives that might be more comfortable for you, see them [here](doc/installation.md). There is a CMake-based building system for several example and test executable files that automatically downloads external dependencies and compiles all executables but it is not needed for the libray's usage.

## Usage

Integrating a function in a specific n-dimensional range is rather simple. You need the following information:
- An *integrator*, a numerical integration method, for which there are [several to choose from](doc/integrators.md).
- An *integrand*, a function to be integrated. It's only parameter has to be a `std::array<F,N>`, where `F` is a floating point number and `N` is the number of dimensions. There are [several ways in which such integrand can be defined](doc/integrands.md).
- A *range*, the integration domain, that is composed on two `std::array<F,N>` marking the limits of the potentially multidimensional integration domain, which [can be defined in different ways](doc/ranges.md).

Example:

```cpp
float sphere(const std::array<float,3>& x) {
    return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<=1.0f?1.0f:0.0f;
}
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<method.integrate(sphere,range)<<"\n";
}
```










