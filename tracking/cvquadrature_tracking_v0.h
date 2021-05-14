#include "region.h"
#include "monte-carlo.h"
#include "integrate-bins-stepper.h"





template<typename Nested, typename Error>
class CVQuadratureTracking {
    
    StepperAdaptive<Nested,Error> stepper;
    unsigned long adaptive_iterations;
    
    template<typename R>
    class Function {
        std::vector<R> regions;
    public:
        template<typename Float>
        typename R::value_type operator()(const std::array<Float,R::dimensions>& x) const {
            for (const auto& r : regions) if (r.range().is_inside(x)) return r.approximation_at(x);
            //Default behavior: extrapolate first region if it is outside all of them
            return regions.front().approximation_at(x);
        }
        template<typename F, typename Float, std::size_t DIM>
        typename R::value_type integral(const F& f, const Range<Float,DIM>& range) const {
            auto sol = regions.front().integral_subrange(range.intersection(regions.front().range()));
            for (auto it = regions.begin()+1; it != regions.end(); ++it)
                sol += (*it).integral_subrange(range.intersection((*it).range()));
            return sol;
        }

        const std::vector<R>& get_regions() const { return regions; }

        Function(std::vector<R>&& rs) : regions(std::forward<std::vector<R>>(rs)) { }
    };

    template<typename R>
    Function<R> function(std::vector<R>&& regions) const { return Function<R>(std::forward<std::vector<R>>(regions)); }
public:
    
    template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {
        //First element of tuple is sum, second element is number of samples
        auto regions = stepper.init(f,range);
        
        for (unsigned long i = 0; i<adaptive_iterations;++i) 
            stepper.step(f,range,regions);


        return Samples<decltype(f(range.min(0)))>();
    }




    template<typename F, typename Float, std::size_t DIM>
    auto operator()(const F& f, const Range<Float,DIM>& range) const {
        auto regions = stepper.init(f,range);
        for (unsigned long i = 0; i<adaptive_iterations;++i) stepper.step(f,range,regions);
        return function(std::move(regions));
    }

    

    CVQuadratureTracking(Nested&& nested, Error&& error, unsigned long ai) :
        stepper(std::forward<Nested>(nested), std::forward<Error>(error)), adaptive_iterations(ai) { }
};
