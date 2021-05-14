#include "region.h"
#include "monte-carlo.h"
#include "integrate-bins-stepper.h"





template<typename Nested, typename Error>
class CVQuadratureTracking {
    
    StepperAdaptive<Nested,Error> stepper;
    unsigned long adaptive_iterations;
    

    template<typename R, typename Result>
    struct Samples {
        R regions;
        Result control; 

        Result sumatory;
        unsigned long counter;
        unsigned long counter_queries;
        Samples(std::vector<R> &rs): 
           regions(std::forward<R>(rs)), sumatory(0),counter(0),counter_queries(0) 
        { }
    };




    /*template<typename R>
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
    Function<R> function(std::vector<R>&& regions) const { return Function<R>(std::forward<std::vector<R>>(regions)); }*/
public:
    
    template<typename F, typename Float>
    auto init(const F& f, const Range<Float,1>& range) const {
        //First element of tuple is sum, second element is number of samples
        using StepperData = decltype(stepper.init(f,range));
        Samples<decltype(f(range.min(0)))>()            
        
        auto regions = stepper.init(f,range);
        
        for (unsigned long i = 0; i<adaptive_iterations;++i) 
            stepper.step(f,range,regions);

        Samples<StepperData, decltype(f(range.min(0)))> samples(regions);
        samples.control = stepper.integral(f, range, regions);
        std::sort(regions.begin(), regions.end(),[]
        )

    }


    template<typename F, typename Float, typename Result>
    void step(const F& f, const Range<Float,1>& range, Samples<Result>& samples, bool verbose = false) const {

        Float maj_sigma_t = samples.control;

        std::exponential_distribution<Float> dis(maj_sigma_t);
        std::uniform_real_distribution<Float> udis01(0,1);

        Float Tr = 1;
        Float Tr_res = 1;
        Float t = range.min(0);

        while(true)
        {

            Float dt = dis(rng);
            t += dt;

            if( t>range.max(0) )
                break;

            
            ++ samples.counter_queries;

            Float Pt = f(t)/(f(t)+fabs(maj_sigma_t-f(t)));
            Float Pn = 1-Pt;


            if( verbose )
                printf("Step: %f [%f, %f] - mu_t(t)= %f / %f: T(t) = %f \n", t, range.min(0), range.max(0), f(t), maj_sigma_t, Tr);


            if( udis01(rng) < Pt)
            {    
                Tr = 0;
                break;
            }
            else
            {
                Tr *= (maj_sigma_t-f(t))/(maj_sigma_t*Pn);
            }

        }


        samples.sumatory += Tr;
        samples.sumatory_res += 
        ++samples.counter;
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
