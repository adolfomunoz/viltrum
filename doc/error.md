# `viltrum` - Error heuristics

There are several error heuristics that can be used for (adaptive integration techniques)[integrators.md]

----

```
error_heuristic_default(<error_metric>)
```

directly applies the error metric. 

----

```
error_heuristic_size(<error_metric>,<size_factor>)
```

This approach calculates the error metric per dimension and adds an additional error term proportional to the size of the range for that dimension, where
- `<size_factor>` is the relative weight of the range size for each dimension. Higher values will favour a more uniform subdivision in all dimensions while lower values will prioritize error estimation for the heuristics. Has a default value of 1.e-5, which worked well in our experiments.

----

There are two error metrics:

```
error_metric_absolute()
```

which applies the standard absolute error metric, and

```
error_metric_relative()
```

which applies a relative error metric.

At some point you might decide that you need a specific error metric for your datatype. If that's the case, you can implement your own error metric by implementing a function that takes two parameters of the same type and returns a floating point number. For that, however, you will need to explore the source code and check how these are implemented.


