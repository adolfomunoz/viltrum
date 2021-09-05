# viltrum

VILTRUM: Varied Integration Layouts for arbiTRary integrals in a Unified Manner - A C++17 header-only library that provides a set of numerical integration routines. This library was generated during research for our paper:

<p align="center">

  <h1 align="center"><a href="https://mcrespo.me/publications/primary-space-cv/">Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations</a></h1>

  <a href="https://mcrespo.me/publications/primary-space-cv/">
    <img src="https://mcrescas.github.io/publications/primary-space-cv/figures/socialMedia.png" alt="Logo" width="100%">
  </a>

  <p align="center">
    ACM Transactions on Graphics - 2021
    <br />
    <a href="https://mcrespo.me/"><strong>Miguel Crespo</strong></a>
    ·
    <a href="http://giga.cps.unizar.es/~ajarabo/"><strong>Adrian Jarabo</strong></a>
    ·
    <a href="http://adolfo-munoz.com/"><strong>Adolfo Muñoz</strong></a>
  </p>

  <p align="center">
    <a href='https://mcrescas.github.io/publications/primary-space-cv/data/crespo2021primary.pdf'>
      <img src='https://img.shields.io/badge/Paper-PDF-red?style=flat-square' alt='Paper PDF'>
    </a>
    <a href='https://mcrescas.github.io/publications/primary-space-cv' style='padding-left: 0.5rem;'>
      <img src='https://img.shields.io/badge/Project-Page-blue?style=flat-square' alt='Project Page'>
    </a>
  </p>
</p>


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
- An *integrand*, a function to be integrated. It's only parameter has to be a `std::array<F,N>`, where `F` is a floating point number and `N` is the number of dimensions. There are [several ways in which such integrand can be defined](doc/integrands.md). The integrand can return any numeric value, or, in general, any data type that supports addition and multiplication by a scalar (tested with [Eigen arrays](https://eigen.tuxfamily.org/dox/group__TutorialArrayClass.html) ).
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

The return type of the `integrate` method will be the same return type of the integrand (`float` in the example above).

### Integrating into bins

It is also possible to obtain not only a single integral value, but to obtain integrals in a *regular n-dimensional grid* of cells or *bins*. For this purpose, this library provides *bin integrators*, that is, integrators that are capable of integrating into bins. For using them, in addition to the integrand and range you need:
- A *bin accesor*, that is, a function that given an *n*-dimensional array (where *n* is the dimensionality of the bin regular structure) of positions (of type `std::size_t`) gives read and write access to the binning structure in that specific position.
- A *resolution*, a *n*-dimensional array of `std::size_t` that gives the resolution of the bin structure for each dimension *n*.

Given these common structure:Bin integrators can be used as follows:

```cpp
float slope(const std::array<float,2>& x) { return (x[1]<x[0])?1.0f:0.0f; }

int main(int argc, char **argv) {
    auto integrator_bins = viltrum::integrator_bins_monte_carlo_uniform(1024);
    auto range = viltrum::range(std::array<float,2>{0,0},std::array<float,2>{1,1});
    float output_array[16];
    //...
}
```

bin integrators can be used as follows (1D example):

```cpp
    auto output_array_access = [&output_array] (const std::array<std::size_t,1>& i) -> float& { return output_array[i[0]]; }; 
    integrator_bins.integrate(output_array_access,std::array<std::size_t,1>{16},slope,range);
    for (float f : output_array) std::cout<<f<<" "; std::cout<<std::endl;
```

and can be expanded to two dimensions as follows:

```cpp
    auto output_array_access_2d = [&output_array] (const std::array<std::size_t,2>& i) -> float& { return output_array[4*i[0]+i[1]]; }; 
    integrator_bins.integrate(output_array_access_2d,std::array<std::size_t,2>{4,4},slope,range);
```

Additionally, the `viltrum` library offers a functional version of the same method:
```cpp
    integrate_bins(integrator_bins,output_array_access,std::array<std::size_t,1>{16},slope,range);
    integrate_bins(integrator_bins,output_array_access_2d,std::array<std::size_t,2>{4,4},slope,range);
```

This functional version is able to deduce the bin accesor as well as the resolution for `std::vector` data types, in which there is only one parameter defining the binning structure:
```cpp
    std::vector<float> output_vector(16);
    integrate_bins(integrator_bins, output_vector, slope, range);
    for (float f : output_vector) std::cout<<f<<" ";
    std::cout<<std::endl;
    
    std::vector<std::vector<float>> output_matrix(4,std::vector<float>(4));
    integrate_bins(integrator_bins, output_matrix, slope, range);
    std::cout<<std::endl;
    for (const std::vector<float>& row : output_matrix) {
        for (float f : row) std::cout<<f<<" ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
```

If you want to define your own automatic bin accesors you can take inspiration from [our source code for the case of vectors](quadrature/bins-containers-adaptor.h) but we do not provide documentation for that (you can always use the general version by explicitly indicating bin accesor and resolution independently as illustrated above).

Last, there are [several bin integrators to choose from](doc/binintegrators.md). However, if you want to use exactly the same integrator with the default parameters than we use in [our paper](https://mcrescas.github.io/publications/primary-space-cv/) define the bin integrator as follows:
  
```cpp
unsigned long spp_cv = std::max(1UL,(unsigned long)(spp*(1.0/16.0)));
auto integrator = integrator_optimized_perpixel_adaptive_stratified_control_variates(
    viltrum::nested(viltrum::simpson,viltrum::trapezoidal), // nested rule, order 2 polynomials
    viltrum::error_single_dimension_size(1.e-5), // error heuristic
    spp_cv*bins/(2*std::pow(3, dim-1)), // number of adaptive iterations calculated from the spps
    std::max(1UL,spp-spp_cv) // number of spps for the residual
);
```
where:
- `spp` is the number of samples per pixel.
- `dim` is the number of dimensions to be explored.
- `bins` is the total number of bins (pixels in the case of images) in all dimensions.  


## License

This code is released under the [GPL v3](LICENSE). Additionally, if you are using this code in academic research, we would be grateful if you cited our paper, for which we generated with this source code:

```bibtex
@article{crespo21primary,
    title = {Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations.},
    year = {2021},
    journal = {ACM Transactions on Graphics},
    author = {Crespo, Miguel and Jarabo, Adrian and Mu\~{n}oz, Adolfo},
    volume = {40},
    number = {3},
    issn = {0730-0301},
    url = {https://doi.org/10.1145/3450627},
    doi = {10.1145/3450627},
    issue_date = {July 2021},
    month = jul,
    articleno = {25},
    numpages = {15},
    keywords = {adaptive quadrature, numerical integration, control variates, rendering, piecewise polynomial}
}
```










