# viltrum

VILTRUM: Varied Integration Layouts for arbiTRary integrals in a Unified Manner — a C++17 header-only library that provides a set of numerical integration routines. This library was developed during research for our paper:

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

`viltrum` is a header-only library and has no external dependencies, so installation is simple (no separate compilation required). Clone the repository into a directory accessible from your build system:

```bash
git clone https://github.com/adolfomunoz/viltrum.git
```

Then ensure the cloned `viltrum` directory is on your compiler include path (for example `-I/path/to/viltrum` for `g++`, or add it with `target_include_directories()` in CMake). Include the main header in your sources:

```cpp
#include <viltrum/viltrum.h>
```

Other options (git submodule, `ExternalProject_Add`, or CMake `FetchContent`) are described [here](doc/installation.md). The repository also contains a CMake-based build for examples and tests; this is optional and not required to use the header-only library.

## Usage

Integrating a function over an n-dimensional domain is straightforward. You need:
- An *integrator*: a numerical integration method (see [integrators](doc/integrators.md)).
- An *integrand*: the function to integrate (see [integrands](doc/integrands.md)). The integrand may return any numeric type that supports addition and scalar multiplication (we tested with [Eigen arrays](https://eigen.tuxfamily.org/dox/group__TutorialArrayClass.html)).
- A *range*: the integration domain and coordinate type (see [ranges](doc/ranges.md)). The range dimensionality must match the integrand's parameters and the numeric type must match the integrand's coordinate type.

The following example shows all three in action: a Monte-Carlo integrator with 64 samples and a random seed, a 1D integrand defined as a lambda function that calculates the sine of a number in `float`, and a 1D integration range between 0 and pi. 

```cpp
float sol = viltrum::integrate(
  viltrum::monte_carlo(64),                         // Monte Carlo with 64 samples
  [] (float x) -> float { return std::sin(x); },    // integrand: sin(x)
  viltrum::range(0.0f, 3.14159265f)                 // range: from 0 to pi
);
```

The return type of `viltrum::integrate` is the integrand's return type (`float` in the example). See the full example [here](../main/doc/montecarlo-1d.cc).

You can also compute integrals on a regular n-dimensional grid of cells (bins). `viltrum` provides bin integrators that accumulate integrals per cell. In addition to the integrand and range you provide:
- A *bin accessor*: a callable that, given an `std::array<std::size_t,N>` bin index, returns a reference to the bin value (for reading/writing).
- A *resolution*: an `std::array<std::size_t,N>` specifying the number of bins along each dimension.

The following example shows how to integrate into a basic array of floats, integrating a 2D integrand into a 1D binning structure. For that, we need to define how to access such array of floats (the *bin accessor*):

```cpp
float sol[10]; // store the results for 10 bins

// Bin accessor: returns a reference into the bin array for a given bin index
auto sol_access = [&sol] (const std::array<std::size_t,1>& pos) -> float& { return sol[pos[0]]; };

// The integrand: f(x,y) = x^2 + y^2 (two-dimensional input)
auto integrand = [] (const std::array<float,2>& x) -> float { return x[0]*x[0] + x[1]*x[1]; };

// Integration range: from (0,0) to (1,1)
auto range = viltrum::range(std::array<float,2>{0.0f, 0.0f}, std::array<float,2>{1.0f, 1.0f});

viltrum::integrate(
  viltrum::monte_carlo(8192),
  sol_access,
  std::array<std::size_t,1>{10}, // resolution (10 bins along the single bin dimension)
  integrand,
  range
);
```

The bin value type (the reference returned by the accessor) and the integrand return type must be compatible. The coordinate type of the range must match the integrand's coordinate type.

Some data types, such as `std::vector` have a specific version of the integrate function that automatically establishes the accessor and the resolution:

```cpp
std::vector<float> sol_vec(10);
viltrum::integrate(
  viltrum::monte_carlo(8192),
  sol_vec, // vector binning structure: accessor and resolution are implicit
  integrand,
  range
);
```

See [main/doc/montecarlo-2d.cc](../main/doc/montecarlo-2d.cc) for the complete example.

Finally, `viltrum` can handle integrals of effectively "infinite" (unbounded) dimensionality — for example, infinite series or some path-integral formulations. Not all integrators support infinite-dimensional integrands, but Monte Carlo-based integrators do. To use infinite-dimensional integration:

- choose an infinite range (e.g. `range_infinite`);
- write the integrand to accept a sequence-like iterable (an automatically generated sequence of samples).

An example:

```cpp
// This integrand is an infinite sum using a Russian-Roulette termination rule.
auto integrand = [] (const auto& seq) -> float {
  auto x = seq.begin();
  const float decaying_factor = 0.5f;
  float sum = 0.0f;
  float term = 1.0f;
  while ((*x) < decaying_factor) {
    ++x;
    term *= ((*x) / decaying_factor);
    ++x;
    sum += term;
  }
  return 2.0f * sum;
};

float sol = viltrum::integrate(
  viltrum::monte_carlo(8192),
  integrand,
  viltrum::range_infinite<float>(0.0f, 1.0f) // infinite-dimensional range (primary domain)
);

// The integral of this function is a geometric series with ratio 0.5, so the result should be ~2.0f.
std::cout << "Integral: " << std::setprecision(6) << sol << " should be close to 2.0\n";
```

See the full example at [main/doc/montecarlo-infd.cc](../main/doc/montecarlo-infd.cc).


## License

This code is released under the [GPL v3](LICENSE). This project includes code from [`pcg-cpp`](https://github.com/imneme/pcg-cpp) and [`Xoshiro-cpp`](https://github.com/Reputeless/Xoshiro-cpp), which are licensed under the MIT License. See the respective header files for detail, as well as the [license for `pcg-cpp`](https://github.com/imneme/pcg-cpp/blob/master/LICENSE-MIT.txt) and the [license for `Xoshiro-cpp`](https://github.com/Reputeless/Xoshiro-cpp/blob/master/LICENSE).


If you use this code in academic research, please consider citing our paper:

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










