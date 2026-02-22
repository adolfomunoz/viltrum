# `viltrum` - Finite dimensionality integrand definition

Defining integrands of finite dimensionality in `viltrum` for their integration with any numerical method is simple. In general, you need to define a function that has the form:

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
    auto method = viltrum::monte_carlo(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<viltrum::integrate(method,sphere,range)<<std::endl;
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
    auto method = viltrum::monte_carlo(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<viltrum::integrate(method,Sphere(),range)<<std::endl;
}
```

- Through a lambda expression:

```cpp
#include <viltrum/viltrum.h>
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::monte_carlo(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<viltrum::integrate(method, [] (const std::array<float,3>& x) { return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1.0f?1.0f:0.0f; },range)<<std::endl;
}
```

## Wrapping normal functions

It is true that defining functions having `std::array`s as parameters is rather unusual and unintuitive. More often you will define a function as having multiple independent parameters. Those are automatically wrapped as `viltrum` integrands, anyway, up to six parameters.

Lets see how this applies with similar examples than above:

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
    std::cout<<"Sphere volume = "<<viltrum::integrate(method,sphere_parameters,range)<<std::endl;
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
    std::cout<<"Sphere volume = "<<viltrum::integrate(method, SphereParameters(),range)<<std::endl;
}
```

- With a lambda expression:
```cpp
#include <viltrum/viltrum.h>
int main() {
    unsigned long samples = 1024;
    auto method = viltrum::integrator_monte_carlo_uniform(samples);
    auto range = viltrum::range(std::array<float,3>{-1.0f,-1.0f,-1.0f},std::array<float,3>{1.0f,1.0f,1.0f});
    std::cout<<"Sphere volume = "<<viltrum::integrate(method,[] (float x, float y, float z) { return (x*x + y*y + z*z)<=1.0f?1.0f:0.0f; },range)<<std::endl;
}
```

The code illustrated above can be tested and compiled in a [source code example](../main/doc/integrands.cc)


# `viltrum` - Infinite dimensionality integrand definition

Defining integrands of "infinite" dimensionality requires passing as parameter an infininte sequence. This infinite sequence is of unknown datatype, so it must be a template. The data type will depend on the numerical integration technique used and cannot be anticipated. Not all integration methods are able to handle integrands of infinite dimensionality.

The most direct way of defining such an integrand is through a class, with a templated `operator()`, as follows:

```cpp
class TaylorExp {
    float rr_prob;
public:
    TaylorExp(float rr_prob = 0.75f) : rr_prob(rr_prob) {}
    template<typename Seq>
    float operator()(const Seq& seq) const {
        float sum = 0.0f; float term = 1.0f; unsigned int n = 0; float prob = 1.0f;   
        auto it = seq.begin(); float x = *it; ++it;
        while ((*it) < rr_prob) {
            prob *= rr_prob; //The probability of reaching this term is the probability of having reached the previous one times the probability of continuing, which is rr_prob.
            sum += term/prob; // sum = sum + term / probability of reaching this term
            ++it; ++n;
            term *= x/float(n); //The next term is the previous one times x/n, so we can update it iteratively. We also use it as a probability, so we don't need to divide by n!
        } 
        return sum;
    }
};
int main(int argc, char** argv) {
    auto method = viltrum::monte_carlo(1024);
    auto range = viltrum::range_primary_infinite<float>();
    std::cout<<viltrum::integrate(method,TaylorExp(),range)<<std::endl;
}
```

Note that, as `operator()` is parametrized through the template parameter `Seq`, it can actually have as parameter any data type that has the same functionality. `Seq` is treated as an infinite sequence from which we can obtain an iterator with the method `begin()`. The iterator has the same functions as a regular iterator (of a `std::vector<float>` for instance).

It is not possible to use a regular templatized C++ function as parameter, because for passing the function itself as parameter you would need to instantiate the template, choosing the parameter in advance, which is currently not possible. However, you can also use lambda functions with an `auto` parameter as follows:

```cpp
int main(int argc, char** argv) {
    auto method = viltrum::monte_carlo(1024);
    auto range = viltrum::range_primary_infinite<float>();
    std::cout<<viltrum::integrate(method,[] (const auto& seq) {
        const float rr_prob = 0.75f; 
        float sum = 0.0f; float term = 1.0f; unsigned int n = 0; float prob = 1.0f;   
        auto it = seq.begin(); float x = *it; ++it;
        while ((*it) < rr_prob) {
            prob *= rr_prob; //The probability of reaching this term is the probability of having reached the previous one times the probability of continuing, which is rr_prob.
            sum += term/prob; // sum = sum + term / probability of reaching this term
            ++it; ++n;
            term *= x/float(n); //The next term is the previous one times x/n, so we can update it iteratively. We also use it as a probability, so we don't need to divide by n!
        } 
        return sum;
    },range)<<std::endl;
}
```

The examples above can be found and compiled from [here](../main/doc/integrands-infinite.cc)