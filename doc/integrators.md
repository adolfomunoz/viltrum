# `viltrum` - Integrators

A `viltrum` integrator represents al algorithm, a numerical method that is able to approximate an integral. All integrators can be used with all different `viltrum::integrate` versions. For instance, for without binning (a single output) the call is as follows:

```cpp
auto solution = viltrum::integrate(<integrator>,<integrand>,<range>) 
```

that applies the algorithm defined on the integrator `<integrator>` to the integral of `<integrand>` along the range `<range>`. There are several ways in which (integrands can be defined)[integrands.md] and in which (integration ranges can be defined)[ranges.md].

When using binning (integrating into a set of bins distributed uniformly in a number of dimensions) then the call can be done as

```cpp
viltrum::integrate(<integrator>,<bins>,<binresolution>,<integrand>,<range>)
```

where `<integrator>`, `<integrand>` and `<range>` have the same meaning. `<binresolution>` is a `std::array<std::size_t,N>` where `N` is the number of dimensions of the binning structure, and each component of the array is the resolution of the binning structure. For instance `std::array<std::size_t,2>{1920,1080}` indicates that the binning structure is two dimensional and has a resolution of 1920x1080. `<bins>` is the binning structure, which is accessed as a function with a `std::array<std::size_t,N>` as parameter which indicates the position for all `N` dimensions (indexed from 0 to the resolution for that dimension minus 1).

There are also specific binning versions for `std::vector<N>` or `std::vector<std::vector<N>>` (where `N` represents any which automatically create the binning structure and the resolution:

```cpp
viltrum::integrate(<integrator>,std::vector<N>,<integrand>,<range>)
viltrum::integrate(<integrator>,std::vector<std::vector<N>>,<integrand>,<range>)
```

For each integrator, its parameters are defined when constructing it, and the type and meaning of each of the parameters depends on the specific integrator. 

For convenience, for the following per-integrator examples we have defined:
- `std::vector<float> sol_bins` as a binning structure, to showcase the fact that all integrators work both with integrating into a single value or into a regular binning structure.
- `integrand` as a finite dimensionality integrand and `range` as its corresponding integration range. [See more about integrand definition](integrands.md)
- `integrand_infinite` as an infinite dimensionality integrand and `range_infinite` as its corresponding range. [See more about how to define integrands of infinite dimensionality](integrands.md)

For the specific implementation of the integrands and their corresponding ranges, you can check the [source code that contains all these examples below](../main/doc/integrators.cc).

## Monte Carlo integrator

[Monte Carlo integration](https://en.wikipedia.org/wiki/Monte_Carlo_integration) uses random numbers to approximate a definite integral, randomly selecting the positions of the evaluated samples within the integration range. The more samples, the longer computation time but the more accuracy.

It can be defined in two ways:

```cpp
viltrum::monte_carlo(<rng>,<nsamples>)
```

where:
- `<rng>` is a C++11 random number generator, such as `std::mt19937` (which is the default) or `std::ranlux48`
- `<nsamples>` which is the number of random samples chosen randomly and uniformly within the integration range.

Alternatively, another construction of this integrator selects by default the `std::mt19937_64` C++11 random number generator as follows:

```
monte_carlo(<nsamples>,<seed>)
```

that will use the default random number generator with the seed `<seed>`. The `<seed>` parameter is optional and, if omitted, it is chosen randomly using C++11's default engine.

Those construction options and parameters are illustrated in the following C++ snipplet:

```cpp
std::cout<<viltrum::integrate(viltrum::monte_carlo(std::ranlux48(),100),integrand, range)<<" ";
std::cout<<viltrum::integrate(viltrum::monte_carlo(100,0), integrand_infinite, range_infinite)<<" ";
viltrum::integrate(viltrum::monte_carlo(100), sol_bins, integrand, range);
```

- The first invocation defines a Monte Carlo integrator with 100 samples and the `std::ranlux48` random number generator.
- The second invocation defines a Monte Carlo integrator with 100 samples and a `std::mt19937` random number generator with seed 0.
- The third invocation defines a Monte Carlo integrator with 100 samples and a `std::mt19937` random number generator with random seed, and integrates using Monte Carlo towards a binning structure.

## Per-bin integrator

The *per-bin integrator* stratifies the integration range into the bins of the binning structure, and then applies another integrator to each of the bins.

```
integrator_per_bin(<bin_integrator>)
```
where `<bin_integrator>` represents the integrator that is applied to each of the bins. There is also a parallel version:
```
integrator_per_bin_parallel(<bin_integrator>)
```
that uses C++17's execution policies. In some compilers (g++) that requires linking [Intel's TBB library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html). An example of the use of this integrator is the following:

```cpp
viltrum::integrate(viltrum::integrator_per_bin(viltrum::monte_carlo(20,0)), sol_bins, integrand, range);
viltrum::integrate(viltrum::integrator_per_bin_parallel(viltrum::monte_carlo(20,0)), sol_bins, integrand_infinite, range_infinite);
```
where the first call integrates into a binning structure using a per-bin integrator that applies Monte Carlo with 20 samples and seed 0 to each of the bins, and the second call does the same but parallelizing the exploration of the bins and applying it to an integrand of infinite dimensionality, which can only be done because the inner integrator supports infinite dimensionality. The use of this integrator only makes sense when integrating into a binning structure.

## Newton-Cotes quadrature rules

[Newton-Cotes formulas](https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas) are a group of formulas that estimate the integral by evaluating a function at regularly spaced sample points, approximating the integrand by a polynomial. Higher order rules are theoretically more accurate than low order rules.

It can be constructed as follows:

```
integrator_newton_cotes(<rule>)
```

where `<rule>` is the specific rule to applied, with the following options:
- `trapezoidal`: the trapezoidal rule (order 1).
- `simpson`: Simpson's rule (order 2).
- `boole`: Boole's rule (order 4).
- `steps<N>(<baserule>)`: composite rule, composed of `N` steps of the `<baserule>` rule.

These can be seen in the following C++ code:

```cpp
std::cout<<viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::trapezoidal), integrand, range)<<" ";
std::cout<<viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::simpson), integrand, range);
std::cout<<viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::boole), integrand, range)<<" ";
viltrum::integrate(viltrum::integrator_newton_cotes(viltrum::steps<2>(viltrum::boole)), sol_bins, integrand, range);
```

where
- The first invocation applies the trapezoidal rule.
- The second invocation applies Simpson's rule.
- The third invocation applies Boole's rule.
- The fourth invocation applies a composite rule of two steps of Boole's rule and integrates it into bins.

None of this integrators works for integrands of infinite dimensionality.

## Adaptive nested Newton-Cotes rules (tolerance-based)

An interesting approach for adaptive integration is the use of nested quadrature rules, which consist of two rules of different order but sharing the same sample points from which the higher order rule estimates the integral while the difference between both rules estimates the error. 

Nested rules enable adaptive integration strategies. By defining a tolerance parameter, integration ranges for which the error exceeds the tolerance are recursively subdivided in half until all subranges yield an error below the tolerance.

The integrator corresponding to this adaptive integration strategy is constructed as follows:

```
integrator_adaptive_tolerance(<nested>,<error>,<tolerance>)
```

where:
- `<nested>` represents a nested quadrature rule, constructed as `nested(<rulehigh>,<rulelow>)` where `<rulehigh>` and `<rulelow>` are the high and low order quadrature rules.
- `<error>` is a [error heuristic](error.md), by default (if omitted) being an absolute error metric.
- `<tolerance>` is a real number that defines the tolerance (should be positive). Lower values yield more accuracy but higher computation time. The specific meaning of the tolerance depends on the error metric. By default is 0.001

This is illustrated in the following C++ code:

```cpp    
std::cout<<viltrum::integrate(
    viltrum::integrator_adaptive_tolerance(
        viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),1.e-3), 
    integrand, range)<<" "; 

viltrum::integrate(
    viltrum::integrator_adaptive_tolerance(
        viltrum::nested(viltrum::boole,viltrum::simpson),
        viltrum::error_heuristic_size(viltrum::error_metric_relative(),1.e-5),1.e-3), 
    sol_bins, integrand, range);
```

where:
- The first line creates an adaptive integrator with tolerance 0.001, with a nested Simpson-Trapezoidal rule and an absolut error metric.
- The second line creates an adaptive integrator with tolerance 0.001 and a nested Boole-Simpson rule, with an heuristic that weights a relative error metric with the size of the region to subdivide with weight 0.00001. Furthermore, the second integrates into bins.

This integrator does not work for integrands of infinite dimensionality.


## Adaptive nested Newton-Cotes rules (iteration-based)

The main problem of a tolerance-based adaptive approach is that the tolerance parameter is heavily linked with the error metric, and it is impossible to anticipate calculation time from such combination. Another option is to order all subranges into a heap according to their estimated error, and keep subdividing the top of the heap and reintroducing into the heap the subranges. This enables a finer control over the computational budget at an additional cost of heap removal and insertion, which is logarithmic with respect to the number of iterations and pays of if this cost is negligible compared with the cost of evaluating the integrand. Each iteration has the same theoretical cost (plus heap insertion/removal) because the number of samples (points in which the integral is evaluated) is proportional to the number of iterations. 

The integrator that applies this strategy is:

```
integrator_adaptive_iterations(<nested>,<error>,<iterations>)
```

where:
- `<nested>` represents a nested quadrature rule, constructed as `nested(<rulehigh>,<rulelow>)` where `<rulehigh>` and `<rulelow>` are the high and low order quadrature rules.
- `<error>` is a (error heuristic)[error.md], by default (if omitted) being an absolute error metric.
- `<iterations>` is a positive integral number that defines the number of iterations. The computational cost is proportional to this number (plus a small extra for heap insertion/removal). 

The following C++ code illustrate this integrator:

```cpp
std::cout<<viltrum::integrate(
    viltrum::integrator_adaptive_iterations(
        viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),32), integrand, range)<<" "; 
viltrum::integrate(
    viltrum::integrator_adaptive_iterations_parallel(
        viltrum::nested(viltrum::boole,viltrum::simpson),
        viltrum::error_heuristic_size(viltrum::error_metric_relative(),1.e-5),32), sol_bins, integrand, range);
```

where:
- The first line creates an adaptive integrator with a nested Simpson-Trapezoidal rule, an absolute error metric and 32 iterations.
- The second line creates an adaptive integrator with a nested Boole-Simpson rule, a error heuristic that weights region size and a relative error metric, and 32 iterations. It is applied for a number of bins, parallelizing the bin exploration. This might require linking with [Intel's TBB library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html). 

This integrator does not work for integrands of infinite dimensionality.

## Fubini integrator

Several of the numerical integrators presented above are unable to work for integrands of infinite dimensionality. Also, there are some numerical integration problems for which specific integrators work better for some dimensions while others work better for other dimensions. To tackle both, we have defined an integrator, based on [Fubini's theorem](https://en.wikipedia.org/wiki/Fubini%27s_theorem) that combines two integrators into a single one, splitting the dimensions. 

The integrator is defined as:

```
fubini<DIM>(lowdim_integrator,highdim_integrator)
```
where
- `DIM` is the number of dimensions covered by the `lowdim_integrator`.
- `lowdim_integrator` is the integrator covers the first `DIM` dimensions.
- `highdim_integrator` is the integrator that covers all dimensions from `DIM`+1 until the last one (even infinite).

An example of the use of this integrator is the following:

```cpp
std::cout<<viltrum::integrate(
    viltrum::integrator_fubini<1>(
        viltrum::integrator_adaptive_iterations_parallel(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),32),
        viltrum::monte_carlo(32)
    ), integrand, range)<<" ";
std::cout<<viltrum::integrate(
    viltrum::integrator_fubini<1>(
        viltrum::integrator_adaptive_iterations_parallel(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_heuristic_default(viltrum::error_metric_absolute()),32),
        viltrum::monte_carlo(32)
    ), integrand_infinite, range_infinite)<<"\n";
```

Both integrators are the same, combined through `integrator_fubini`, were the first dimension is an adaptive newton cotes and the rest of the dimensions is Monte-Carlo. The upper shows the use on an integrand of finite dimensionality and the second of infinite dimensionality.


## Adaptive nested Newton-Cotes control variates with Monte Carlo integration of the residual

This strategy is the base of our paper [**Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations**](https://mcrescas.github.io/publications/primary-space-cv/), and it preserves the best of both strategies: the low frequency regions are better recovered using adaptive Newton-Cotes for a number of iterations and high frequency details are better recovered using Monte-Carlo (of the residual with respect to the Newton-Cotes approximation). 

This integrator follows a sampling strategy for the residual in which each individual region, represented by a single multivariate polynomial in the piecewise control variate, is selected randomly and uniformly, and then again uniformly inside each region. This means that smaller regions (with allegedly larger error estimations) will get more sample density than larger regions. 

This integrator is constructed as follows:

```
integrator_crespo21(<iterations>,<spp>,<seed>)
```
where
- `<iterations>` is a positive integral number that defines the number of iterations for constructing the piecewise-polynomial control variate. The computational cost of this construction is proportional to this number (plus a small extra for heap insertion/removal). 
- `<spp>` which is the number of random samples chosen within the integration range following the subdivisions of the piecewise polynomials.
- `<seed>` is the seed of the random number generator. If omitted it is calculated randomly.

There is a version for integrands of infinite dimensionality
```
integrator_crespo21_infinite<DIM>(<iterations>,<iteration_samples>,<spp>,<seed>)
```
where all items represent the same thingy, `<DIM>` represents the number of dimensions covered by the polynomial, and the new one, `<iteration_samples>`, represents the number of Monte Carlo samples used for estimating the higher dimensions of each sample required to calculate the polynomial.

This is illustrated in the following code:

```cpp
std::cout<<viltrum::integrate(
    viltrum::integrator_crespo2021(16,64), integrand, range)<<" ";
std::cout<<viltrum::integrate(
    viltrum::integrator_crespo2021_infinite<4>(16,4,64), integrand_infinite, range_infinite)<<"\n";
```

The code illustrated in this page can be tested and compiled in a [source code example](../main/doc/integrators.cc)


