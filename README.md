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
    <strong>Adrian Jarabo</strong>
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
- An *integrand*, a function to be integrated. There are [several ways in which such integrand can be defined](doc/integrands.md). The integrand can return any numeric value, or, in general, any data type that supports addition and multiplication by a scalar (tested with [Eigen arrays](https://eigen.tuxfamily.org/dox/group__TutorialArrayClass.html) ).
- A *range*, the integration domain, that marks the limits of the integration domain, which [can be defined in different ways](doc/ranges.md). The dimensionality of the integration range must match the dimensionality of the integrand, and the numeric data type of the range must match the data type of the parameters.

The following example shows all three in action: a Monte-Carlo integrator with 64 samples and a random seed, a 1D integrand defined as a lambda function that calculates the sine of a number in `float`, and a 1D integration range between 0 and pi. 

```cpp
float sol = viltrum::integrate(
    viltrum::monte_carlo(64),  //Numerical integration techinque: Monte-Carlo with 64 samples
    [] (float x) -> float { return std::sin(x);}, //Function to integrate: sin(x)
    viltrum::range(0.0f,3.14159265f) //Integration range: from 0 to pi
);
```

The return type of the `integrate` function will be the same return type of the integrand (`float` in the example above). The full simple code for this example is [here](main/doc/montecarlo-1d.cc).

Alternatively, it is also possible to obtain not only a single integral value, but to obtain integrals in a *regular n-dimensional grid* of cells or *bins*. For this purpose, this library provides *bin integrators*, that is, integrators that are capable of integrating into bins. For using them, in addition to the integrand and range you need:
- A *bin accesor*, that is, a function that given an *n*-dimensional array (where *n* is the dimensionality of the bin regular structure) of positions (indices of type `std::size_t`) gives read and write access to the binning structure in that specific position.
- A *resolution*, a *n*-dimensional array of `std::size_t` that gives the resolution of the bin structure for each dimension *n*.

The following example shows how to integrate into a basic array of floats, integrating a 2D integrand into a 1D binning structure. For that, we need to define how to access such array of floats (the *bin accessor*):

```cpp
float sol[10]; //We will compute the integral in 10 bins, so we need an array of 10 floats to store the results
//Bin accesor: access to the binning structure (sol) through the position (pos, the index) in the range. 
auto sol_access = [&sol] (const std::array<std::size_t,1>& pos) -> float& { return sol[pos[0]]; }; 
//  The function to integrate is f(x,y) = x^2 + y^2. The parameter is an std::array with the two dimenions but it could also be a bidimensional function with two parameters
auto integrand = [] (const std::array<float,2>& x) -> float { return x[0]*x[0] + x[1]*x[1]; };
//Integration range: from (0,0) to (1,1). Could also be viltrum::range(0.0f,1.0f,0.0f,1.0f)
auto range = viltrum::range(std::array<float,2>{0.0f,0.0f}, std::array<float,2>{1.0f,1.0f}); 
viltrum::integrate(
    viltrum::monte_carlo(8192),  //Numerical integration techinque: Monte-Carlo with 8192 samples (distributed along the whole integration domain)
    sol_access,
    std::array<std::size_t,1>{10}, //Number of bins in each dimension (8 bins in this case) integrand. Dimensionality should be the same than the parameter of the access above
    integrand,
    range
);
```

The data type of the binning structure (the reference returned by the bin accesor) and the data type returned by the integrand must be the same (or at least, compatible). Also, as above, the data type of each element of the range (`float` in the example above) must be of the same type of every element of the parameter set of the function. The full code of the example above is [here](main/doc/montecarlo-2d.cc).

Last, it is also possible to deal with integrals of "infinite" (unbounded) dimensionality. Examples of such are an infinite series or the path integral for rendering. Not all integrators are able to deal with integrals of unbounded dimensionality (Monte-Carlo, deals with it perfectly, though). Integrals of unbounded dimensionality are dealt with easily: first, the integration range must be of infinite dimensionality (`range_infinite`). Then, the integrand must have as parameter an "automatic" data type, which is an iterable infinite sequence of numbers (the paramters of the integrand). An example of this is the following:

```cpp
// This integrand is an infinite sum of decaying products, with decaying factor 0.5.
auto integrand = [] (const auto& seq) -> float {
    auto x = seq.begin(); //The sequence is an infinite sequence of numbers within the integration range, which is generated on the fly by the integrator. We can iterate over it until we want.
    const float decaying_factor = 0.5f;
    float sum = 0.0f;
    float term = 1.0f;  
    // For each iteration we take two samples:
    // - the first is Russian Roulette for the component with the decaying factor.
    // - the second is the value of the component with the decaying factor, which is multiplied by the term and added to the sum.       
    while ((*x) < decaying_factor) { 
        ++x; 
        term *= ((*x)/decaying_factor); 
        ++x; 
        sum += term; 
    }
    return 2.0f*sum;
};
float sol = viltrum::integrate(
    viltrum::monte_carlo(8192),  //Numerical integration techinque: Monte-Carlo with 8192 samples
    integrand,
    viltrum::range_infinite<float>(0.0f,1.0f) //Range of infinite dimensionality, all dimensions between 0 and 1
);
// The integral of this function is a geometric series with ratio 0.5, so the result should be 2.0f.
std::cout<<"Integral: "<<std::setprecision(6)<<sol<<" should be close to 2.0\n";
```

You can find the full example [here](main/doc/montecarlo-infd.cc).


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










