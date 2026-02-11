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

where `<integrator>`, `<integrand>` and `<range>` have the same meaning. `<binresolution>` is a `std::array<std::size_t,N>` where `N` is the number of dimensions of the binning structure, and each component of the array is the resolution of the binning structure. For instance `std::array<std::size_t,2>{1920,1080}` indicates that the binning structure is two dimensional and has a resolution of 1920x1080. `<bins>` is the binning structure, which is accessed as a function with a `std::array<std::size_t,N>` as parameter.

For each integrator, its parameters are defined when constructing it, and the type and meaning of each of the parameters depends on the specific integrator.

## Monte Carlo integrator

[Monte Carlo integration](https://en.wikipedia.org/wiki/Monte_Carlo_integration) uses random numbers to approximate a definite integral, randomly selecting the positions of the evaluated samples within the integration range. The more samples, the longer computation time but the more accuracy.

It can be defined in two ways:

```
integrator_monte_carlo_uniform(<rng>,<nsamples>)
```

where:
- `<rng>` is a C++11 random number generator, such as `std::mt19937_64` (which is the default) or `std::ranlux48`
- `<nsamples>` which is the number of random samples chosen randomly and uniformly within the integration range.

Alternatively, another construction of this integrator selects by default the `std::mt19937_64` C++11 random number generator as follows:

```
integrator_monte_carlo_uniform(<nsamples>,<seed>)
```

that will use the default random number generator with the seed `<seed>`. The `<seed>` parameter is optional and, if omitted, it is chosen randomly using C++11's default engine.

Those construction options and parameters are illustrated in the following C++ snipplet:

```cpp
std::cout<<viltrum::integrator_monte_carlo_uniform(std::ranlux48(),100).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_monte_carlo_uniform(100,0).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_monte_carlo_uniform(100).integrate(function,range)<<"\n";
```

- The first invocation defines a Monte Carlo integrator with 100 samples and the `std::ranlux48` random number generator.
- The second invocation defines a Monte Carlo integrator with 100 samples and a `std::mt19937_64` random number generator with seed 0.
- The third invocation defines a Monte Carlo integrator with 100 samples and a `std::mt19937_64` random number generator with random seed.

## Newton-Cotes quadrature rules

[Newton-Cotes formulas](https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas) are a group of formulas that estimate the integral by evaluating a function at regularly spaced sample points, approximating the integrand by a polynomial. Higher order rules are theoretically more accurate than low order rules.

It can be constructed as follows:

```
integrator_quadrature(<rule>)
```

where `<rule>` is the specific rule to applied, with the following options:
- `trapezoidal`: the trapezoidal rule (order 1).
- `simpson`: Simpson's rule (order 2).
- `boole`: Boole's rule (order 4).
- `steps<N>(<baserule>)`: composite rule, composed of `N` steps of the `<baserule>` rule.

These can be seen in the following C++ code:

```cpp
std::cout<<viltrum::integrator_quadrature(viltrum::trapezoidal).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_quadrature(viltrum::simpson).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_quadrature(viltrum::boole).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_quadrature(viltrum::steps<2>(viltrum::boole)).integrate(function,range)<<"\n";
```

where
- The first invocation applies the trapezoidal rule.
- The second invocation applies Simpson's rule.
- The third invocation applies Boole's rule.
- The fourth invocation applies a composite rule of two steps of Boole's rule.

## Adaptive nested Newton-Cotes rules (tolerance-based)

An interesting approach for adaptive integration is the use of nested quadrature rules, which consist of two rules of different order but sharing the same sample points from which the higher order rule estimates the integral while the difference between both rules estimates the error. 

Nested rules enable adaptive integration strategies. By defining a tolerance parameter, integration ranges for which the error exceeds the tolerance are recursively subdivided in half until all subranges yield an error below the tolerance.

The integrator corresponding to this adaptive integration strategy is constructed as follows:

```
integrator_adaptive_tolerance(<nested>,<error>,<tolerance>)
```

where:
- `<nested>` represents a nested quadrature rule, constructed as `nested(<rulehigh>,<rulelow>)` where `<rulehigh>` and `<rulelow>` are the high and low order quadrature rules.
- `<error>` is a (error metric)[error.md], by default (if omitted) being an absolute error metric separated by dimension (`error_absolute_single_dimension()`).
- `<tolerance>` is a real number that defines the tolerance (should be positive). Lower values yield more accuracy but higher computation time. The specific meaning of the tolerance depends on the error metric. By default is 0.001

This is illustrated in the following C++ code:

```cpp    
std::cout<<viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),1.e-3).integrate(function,range)<<" ";   
std::cout<<viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::boole,viltrum::simpson),1.e-3).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_tolerance(viltrum::nested(viltrum::steps<2>(viltrum::boole),viltrum::boole)).integrate(function,range)<<"\n";
```

where:
- The first line creates an adaptive integrator with tolerance 0.001, with a nested Simpson-Trapezoidal rule and a relative error metric per dimension.
- The second line creates an adaptive integrator with tolerance 0.001 and a nested Boole-Simpson rule, with the default absolute error metric per dimension.
- The third line creates an adaptive integrator with a nested 2 steps Boole-1 step Boole rule, with the default 0.001 tolerance and the default absolute error metric per dimension. 

## Adaptive nested Newton-Cotes rules (iteration-based)

The main problem of a tolerance-based adaptive approach is that the tolerance parameter is heavily linked with the error metric, and it is impossible to anticipate calculation time from such combination. Another option is to order all subranges into a heap according to their estimated error, and keep subdividing the top of the heap and reintroducing into the heap the subranges. This enables a finer control over the computational budget at an additional cost of heap removal and insertion, which is logarithmic with respect to the number of iterations and pays of if this cost is negligible compared with the cost of evaluating the integrand. Each iteration has the same theoretical cost (plus heap insertion/removal) because the number of samples (points in which the integral is evaluated) is proportional to the number of iterations. 

The integrator that applies this strategy is:

```
integrator_adaptive_iterations(<nested>,<error>,<iterations>)
```

where:
- `<nested>` represents a nested quadrature rule, constructed as `nested(<rulehigh>,<rulelow>)` where `<rulehigh>` and `<rulelow>` are the high and low order quadrature rules.
- `<error>` is a (error metric)[error.md], by default (if omitted) being an absolute error metric separated by dimension (`error_absolute_single_dimension()`).
- `<iterations>` is a positive integral number that defines the number of iterations. The computational cost is proportional to this number (plus a small extra for heap insertion/removal). 

The following C++ code illustrate this integrator:

```cpp
std::cout<<viltrum::integrator_adaptive_iterations(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),10).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_iterations(viltrum::nested(viltrum::boole,viltrum::simpson),10).integrate(function,range)<<"\n";
```

where:
- The first line creates an adaptive integrator with a nested Simpson-Trapezoidal rule, a relative error metric per dimension and 10 iterations.
- The second line creates an adaptive integrator with a nested Boole-Simpson rule, the default error metric and 10 iterations.


<!--
## Adaptive nested Newton-Cotes control variates with Monte Carlo integration of the residual

This strategy is the base of our paper [**Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations**](https://mcrescas.github.io/publications/primary-space-cv/), and it preserves the best of both strategies: the low frequency regions are better recovered using adaptive Newton-Cotes for a number of iterations and high frequency details are better recovered using Monte-Carlo (of the residual with respect to the Newton-Cotes approximation). 

This integrator follows a sampling strategy for the residual in which each individual region, represented by a single multivariate polynomial in the piecewise control variate, is selected randomly and uniformly, and then again uniformly inside each region. This means that smaller regions (with allegedly larger error estimations) will get more sample density than larger regions. 

This integrator is constructed as follows:

```
integrator_adaptive_control_variates(<nested>,<error>,<iterations>,<rng>,<nsamples>)
```

where
- `<nested>` represents a nested quadrature rule, constructed as `nested(<rulehigh>,<rulelow>)` where `<rulehigh>` and `<rulelow>` are the high and low order quadrature rules.
- `<error>` is a (error metric)[error.md], by default (if omitted) being an absolute error metric separated by dimension (`error_absolute_single_dimension()`).
- `<iterations>` is a positive integral number that defines the number of iterations for constructing the piecewise-polynomial control variate. The computational cost of this construction is proportional to this number (plus a small extra for heap insertion/removal). 
- `<rng>` is a C++11 random number generator, such as `std::mt19937_64` (which is the default) or `std::ranlux48`
- `<nsamples>` which is the number of random samples chosen randomly and uniformly within the integration range.

Alternativelly, instead of indicating a random number generator it is possible to use the default `std::mt19937_64` generator by just writing the `<seed>`, as follows:

```
integrator_adaptive_control_variates(<nested>,<error>,<iterations>,<nsamples>,<seed>)
```

where `<seed>` is the seed that, if omitted, is calculated randomly.

These constructions are illustrated in the following code:

```cpp
std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),4,std::ranlux48(),80).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,std::ranlux48(),80).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),viltrum::error_relative_single_dimension(),4,80,0).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,std::ranlux48(),80,0).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,80,0).integrate(function,range)<<" ";
std::cout<<viltrum::integrator_adaptive_control_variates(viltrum::nested(viltrum::simpson,viltrum::trapezoidal),4,80).integrate(function,range)<<"\n";
```

where
- The first call constructs an adaptive control variate using the Simpson-Trapezoidal nested rule, with a relative error metric and four adaptive iterations, and the residual uses 80 Monte Carlo samples using the `std::ranlux48()` random number generator.
- The second call constructs an adaptive control variate using the Simpson-Trapezoidal nested rule, with the default absolute error metric and four adaptive iterations, and the residual uses 80 Monte Carlo samples using the `std::ranlux48()` random number generator.
- The third call constructs an adaptive control variate using the Simpson-Trapezoidal nested rule, with a relative error metric and four adaptive iterations, and the residual uses 80 Monte Carlo samples using the default random number generator with seed 0.
- The fourth call constructs an adaptive control variate using the Simpson-Trapezoidal nested rule, with the default absolute error metric and four adaptive iterations, and the residual uses 80 Monte Carlo samples using default random number generator width seed 0.
- The last call constructs an adaptive control variate using the Simpson-Trapezoidal nested rule, with the default absolute error metric and four adaptive iterations, and the residual uses 80 Monte Carlo samples using default random number generator width random seed.




The code illustrated in this page can be tested and compiled in a [source code example](../main/doc/integrators.cc)

-->

