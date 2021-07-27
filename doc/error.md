# `viltrum` - Error metrics

There are several error metrics that can be used for (adaptive integration techniques)[integrators.md]

----

```
error_absolute_single_dimension()
```

calculates a standard absolute error for each dimension. 

----

```
error_error_single_dimension(<offset>)
```

calculates a standard relative error (absolute error in the numerator with respect to value in the denominator) for each dimension, where
- `<offset>` is a floating point number, added to the value in the denominator in order to avoid divisions by zero. The default value is 1.e-6

----

```
error_absolute_single_dimension_size(<size_factor>)
```

This approach calculates the standard absolute error per dimension and adds an additional error term proportional to the size of the range for that dimension, where
- `<size_factor>` is the relative weight of the range size for each dimension. Higher values will favour a more uniform subdivision in all dimensions while lower values will prioritize error estimation for the heuristics. Has a default value of 1.e-5, which worked well in our experiments.

This error metric has the advantage of using a reasonable absolute error metric but prevent algorithms for getting stagnant due to an incorrect error estimation (because eventually the size will matter for the error estimation). This is the default error metric for our paper (Primary-Space Adaptive Control Variates using Piecewise-Polynomial Approximations)[https://mcrescas.github.io/publications/primary-space-cv/].

----

```
error_relative_single_dimension_size(<size_factor>,<offset>)
```

which is similar to the previous one (includes a term for the size of the range into the error estimation) but also uses a relative error metric, where
- `<size_factor>` is the relative weight of the range size for each dimension. Higher values will favour a more uniform subdivision in all dimensions while lower values will prioritize error estimation for the heuristics. Has a default value of 1.e-5, which worked well in our experiments.
- `<offset>` is a floating point number, added to the value in the denominator in order to avoid divisions by zero. The default value is 1.e-6

