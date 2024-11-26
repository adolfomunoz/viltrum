#pragma once
#include "norm.h"



namespace viltrum {
    class cv_fixed_weight {
    double alpha; //This should go from zero (full importance sampling) to 1 (full control variates)
public:
    cv_fixed_weight(double a = 1) : alpha(a) {}

    template<typename Sample>
    class Accumulator {
    private:
        Sample sum; std::size_t size; double alpha;
        Accumulator(double alpha, const Sample& ini = Sample(0)) : 
            sum(ini), size(0), alpha(alpha) {}
    public:
        void push(const Sample& function_sample, const Sample& approximation_sample) {
            sum += function_sample - alpha*approximation_sample;
            ++size;
        }

        Sample integral(const Sample& approximation) const {
            if (size == 0) return approximation;
            else return sum/double(size) + alpha*approximation;
        }
        friend class cv_fixed_weight;
    };

    template<typename Sample>
    Accumulator<Sample> accumulator(const Sample& ini = Sample(0)) const {
        return Accumulator(alpha,ini);
    }
};

template<typename Norm = NormDefault>
class cv_optimize_weight {
    Norm norm;
public:
    cv_optimize_weight(const Norm& n = Norm()) : norm(n) {}
    template<typename Sample>
    class Accumulator {
    private:
        Norm norm;
        Sample sum_f, sum_app; std::size_t size;
        //These are for online calculation of variance and covariance
        double k_f, k_app, e_f, e_ap, e_ap2, e_fap;
        Accumulator(const Norm& n, const Sample& ini = Sample(0)) : 
            norm(n),sum_f(ini),sum_app(ini),size(0),
            k_f(0), k_app(0), e_f(0), e_ap(0), e_ap2(0), e_fap(0) {}
    public:
        void push(const Sample& function_sample, const Sample& approximation_sample) {
            //For online variance and covariance (and alpha)
            if (size==0) {
                k_f = norm(function_sample); k_app = norm(approximation_sample);
            }

            e_f += (norm(function_sample) - k_f);
            e_ap += (norm(approximation_sample) - k_app);
            e_ap2 += (norm(approximation_sample) - k_app)*(norm(approximation_sample) - k_app);
            e_fap += (norm(function_sample) - k_f)*(norm(approximation_sample) - k_app);

            sum_f += function_sample; sum_app += approximation_sample; 
            ++size;
        }

        double covariance() const {
            return (e_fap - (e_f*e_ap)/double(size))/double(size-1);
        } 

        double variance() const {
            return (e_ap2 - (e_ap*e_ap)/double(size))/double(size-1);
        } 

        double alpha() const {
            if (size < 2) return 1;
            auto c = std::max(0.0,covariance());
            if (c<=0.0) return 0.0;
            auto v = std::max(c,variance());
            return c/v;
        } 

        Sample integral(const Sample& approximation) const {
            if (size < 2) return approximation;
            else {
                auto a = alpha();
    //            std::cerr<<"Alpha = "<<a<<std::endl;
                return (sum_f - a*sum_app)/double(size) + a*approximation;
            }
        }
        friend class cv_optimize_weight;
    };

    template<typename Sample>
    Accumulator<Sample> accumulator(const Sample& ini = Sample(0)) const {
        return Accumulator<Sample>(norm,ini);
    }
};
} 