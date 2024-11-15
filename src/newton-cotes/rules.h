#pragma once
#include <array>
#include <vector>
#include <list>
#include "../control-variates/norm.h"

namespace viltrum {

struct Trapezoidal {	
	static constexpr std::size_t samples = 2;
	
	template<typename T>
	auto operator()(const std::array<T,samples>& p) const -> T {
		return (p[0]+p[1])/2.0;
	}

    template<typename T>
    constexpr T coefficient(const std::array<T,samples>& p, std::size_t i) const {
        switch (i) {
            case 1:  return p[1]-p[0];
            default: return p[0];
        }
    }


    template<typename T>
    constexpr std::array<T, samples> coefficients(const std::array<T,samples>& p) const {
        return std::array<T, samples>{
                  p[0],
                  p[1]-p[0]
        };
    }

    template<typename Float, typename T>	
	constexpr T at(Float t, const std::array<T,samples>& p) const {
        auto c = coefficients(p);
		return c[1]*t + c[0];
	}
	
	template<typename Float, typename T>
	constexpr T subrange(Float a, Float b, const std::array<T,samples>& p) const {
        auto c = coefficients(p);
		return (c[1]*b/2.0 + c[0])*b - (c[1]*a/2.0 + c[0])*a;
	}

/*
	template<typename Float, typename T>	
	auto at(Float t, const std::array<T,samples>& p) const -> T {
		return p[0]*(1-t)+p[1]*t;
	}
	
	template<typename Float, typename T>
	auto subrange(Float a, Float b, const std::array<T,samples>& p) const -> T {
		return b*(p[0] + (p[1] - p[0])*b/2.0) - a*(p[0] + (p[1] - p[0])*a/2.0);  
	}
*/
} trapezoidal;

struct Simpson {
	static constexpr std::size_t samples = 3;
		
	template<typename T>
	T operator()(const std::array<T,samples>& p) const {
		return (p[0]+4.0*p[1]+p[2])/6.0;
	}

    template<typename T>
    constexpr T coefficient(const std::array<T,samples>& p, std::size_t i) const {
        switch (i) {
            case 1: return -3*p[0]+4*p[1]-p[2];
            case 2: return 2*p[0]-4*p[1]+2*p[2];
            default: return p[0];
        }
    }

    template<typename T>
    constexpr std::array<T, samples> coefficients(const std::array<T,samples>& p) const {
        return std::array<T, samples>{
                  p[0],
                  -3*p[0]+4*p[1]-p[2],
                  2*p[0]-4*p[1]+2*p[2]
        };
    }

	template<typename Float, typename T>	
	constexpr T at(Float t, const std::array<T,samples>& p) const {
	    auto c = coefficients(p);
		return (c[2]*t + c[1])*t + c[0];
	}

	static constexpr bool isnan(double d) { return std::isnan(d); }
	static constexpr bool isnan(float d) { return std::isnan(d); }
	template<typename V>
	static constexpr bool isnan(const V& v) { return isnan(v[0]); }

	template<typename Float, typename T>	
	constexpr T subrange(Float a, Float b, const std::array<T,samples>& p) const {
	    auto c = coefficients(p);
		return ((c[2]*b/3.0 + c[1]/2.0)*b + c[0])*b - ((c[2]*a/3.0 + c[1]/2.0)*a + c[0])*a;	
	}

	static constexpr double pdf_epsilon = 1.e-5;

private:
	template<typename Float, typename T, typename Norm>
	std::array<Float,samples> pdf_points(Float t0, Float t1, const std::array<T,samples>& p, const Norm& norm = Norm()) const {
		std::array<Float,samples> p_pdf; 
//		for (std::size_t i = 0; i<samples; ++i) p_pdf[i] = norm.sign(p[i]);
		for (std::size_t i = 0; i<samples; ++i) p_pdf[i] = norm(p[i]); //Only possitive side
/*
		Not anymore, we only take care of the possitive side
		if (subrange(t0,t1,p_pdf)<0) { //If the integral is mostly negative, the pdf tries to take care of the negative side
			for (std::size_t i = 0; i<samples; ++i) p_pdf[i] *= -1;
		}
*/
		std::array<Float,samples> cs = coefficients(p_pdf);

		Float ymin = 0;

		if (cs[2] > 0) { //The minimum might in the middle of the polynomial 
			//We look at the division
			Float tmin = -cs[1]/(2.0*cs[2]);
			if ((tmin > 0) && (tmin < 1)) {
				Float ytmin = (cs[2]*tmin + cs[1])*tmin + cs[0];
				if (ytmin < ymin) ymin = ytmin;
			} 
		}

		//We check all the points as well
		//  for (Float p : p_pdf) if (p<ymin) ymin=p;
		// No need, we are on the possitive side only

/*
		//Particular case: if coefficients and points are really really small we increase ymin to avoid
		//  pdfs that are completelly zero
		Float sum_all = 0;
		for (std::size_t i = 0; i<samples; ++i) sum_all += std::abs(p_pdf[i]-ymin);
		if (sum_all<pdf_epsilon) ymin -= 0.01; //It is already uniform sampling so we just make sure that it is not zero but not big
*/
		//We instead always give a margin so there are never probability 0 which would be problematic
		//ymin-=0.01;

		for (std::size_t i = 0; i<samples; ++i) p_pdf[i] -= ymin;

		return p_pdf;	
	}
	
public:
	//This is the same than above but for the pdf so it integrates the absolute norm of the polynomial
	template<typename Float, typename T, typename Norm = NormDefault>
	constexpr Float pdf_integral_subrange(Float t0, Float t1, 
		const std::array<T,samples>& p, const Norm& norm = Norm()) const {
			return subrange(t0,t1,pdf_points(t0,t1,p,norm));
		}

	//Unnormalized pdf...
	template<typename Float, typename T,typename Norm = NormDefault>	
	constexpr Float pdf_unnormalized(Float t, const std::array<T,samples>& p, Float a, Float b, const Norm& norm = Norm()) const {
		auto ps = pdf_points(a,b,p,norm);
		return at(t,ps);
	}

	//Normalized pdf
	template<typename Float, typename T,typename Norm = NormDefault>	
	constexpr Float pdf(Float t, const std::array<T,samples>& p, Float a, Float b, const Norm& norm = Norm()) const {
		auto ps = pdf_points(a,b,p,norm);
		return at(t,ps)/subrange(a,b,ps);
	}
	//I had to create this function, otherwise c++ returns -nans for a cube root of a negative number, even though 
	//it's a cube root
	template<typename Float>
	constexpr static Float cuberoot(Float x) { 
		return (x<0)?(-1*pow(-x,1./3.)):pow(x,1./3.);
	}

	//Name cdf, inv_cdf and sample_normalized as you wish
	//cdf, inv_cdf and sample_normalize do not work unless there is no roots in the order two polynomial
	template<typename Float, typename T>
	constexpr T cdf(Float t, const std::array<T,samples>& p) const {
    	auto c = coefficients(p);
		return ((c[2]*t/3.0 + c[1]/2.0)*t + c[0])*t;
	}

	template<typename Float, typename T, typename Norm>
	std::vector<Float> inv_cdf(Float x, const std::array<T,samples>& points, const Norm& norm = Norm()) const {
		std::vector<Float> solutions;

		auto coeff = coefficients(points);
		auto a = norm.sign(coeff[2])/3.;	//Term multiplying x³ in the cdf
		auto b = norm.sign(coeff[1])/2.;	//Term multiplying x² in the cdf
		auto c = norm.sign(coeff[0]);		//Term multiplying x in the cdf
		auto d = -x;			//-x because we are solving the roots for cdf = x (cdf - x = 0) to compute the inverse

		//We renormalize a, b, c and d for stability and boundary checking
		auto sum = a+b+c+d; a/=sum; b/=sum; c/=sum; d/=sum;

		if (std::abs(a)<1.e-3) { //Second degree equation	
			if (std::abs(b)<1.e-3) { //Degree one equation
				if (std::abs(c)>=1.e-3) //If all zeroes no solution   
					solutions.push_back(-d/c);
			} else {
				auto disc = c*c - 4*b*d;
				if (disc>=0) {
					solutions.push_back((-c + std::sqrt(disc))/(2*b));
					solutions.push_back((-c - std::sqrt(disc))/(2*b));
				}  
			} 
		} else { 
			//In case you are curious, this is Cardano's method for one solution:
			auto p = c/a - pow(b,2.)/(3.*pow(a,2.));
			auto q = 2*pow(b,3.)/(27.*pow(a,3.)) - b*c/(3.*pow(a,2.)) + d/a;
//			std::cerr<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<p<<" "<<q<<std::endl; 
			auto sqr = pow(q/2.,2.) + pow(p/3.,3.);
			if(sqr >= 0.){ 
				//One solution
				solutions.push_back(cuberoot(-q/2. - sqrt(sqr)) + cuberoot(-q/2. + sqrt(sqr))-b/(3*a));
			}
			else{ 
				//3 solutions (Viète's trigonometric expression of the roots)
				for(int i=0; i<3; i++)
					solutions.push_back(2.*sqrt(-p/3.) * cos(1./3. * acos(3.*q/(2.*p) * sqrt(-3./p)) - 2.*M_PI * i/3.)-b/(3.*a));
			}
		}	
		return solutions;
	}

	template<typename Float, typename T, typename Norm = NormDefault>
	constexpr Float sample(Float s, Float t0, Float t1, const std::array<T,samples>& p,
			const Norm& norm = Norm()) const {

		auto ps = pdf_points(t0,t1,p,norm);
		Float x = s*cdf(t1,ps) + (Float(1)-s)*cdf(t0,ps);

		auto res = inv_cdf(x,ps,norm);
		for (Float r : res) {
			//In some corner cases we are getting a value which is extremelly close to the boundaries so we 
			// have added an "epsilon" to the range check
			if ((r<t0) && (std::abs(t0 - r)<1.e-5)) r = t0;
			if ((r>t1) && (std::abs(r - t1)<1.e-5)) r = t1;
			if ((r>=t0) && (r<=t1) && (!std::isnan(r)))	{ 
				return r; 
			}	
		}

		//In weird cases that should not happen we default to uniform sampling
		return s*(t1-t0) + t0;
	}	

} simpson;

struct Boole {
	static constexpr std::size_t samples = 5;
			
	template<typename T>
	auto operator()(const std::array<T,samples>& p) const -> T {
		return (7.0*p[0]+32.0*p[1]+12.0*p[2]+32.0*p[3]+7.0*p[4])/90.0;
	}

    template<typename T>
    constexpr T coefficient(const std::array<T,samples>& p, std::size_t i) const {
        switch (i) {
            case 1: return -25*p[0]/3+ 16*p[1] - 12*p[2] +16*p[3]/3 - p[4];
            case 2: return 70*p[0]/3 -208*p[1]/3 + 76*p[2] -112*p[3]/3 + 22*p[4]/3;
            case 3: return -80*p[0]/3 + 96*p[1] - 128*p[2] + 224*p[3]/3 - 16*p[4];
            case 4: return 32*p[0]/3 - 128*p[1]/3 + 64*p[2] - 128*p[3]/3 + 32*p[4]/3;
            default: return p[0];
        }
    }


    template<typename T>
    constexpr std::array<T, samples> coefficients(const std::array<T,samples>& p) const {
        return std::array<T, samples>{
                  p[0],
                  -25*p[0]/3+ 16*p[1] - 12*p[2] +16*p[3]/3 - p[4],
                  70*p[0]/3 -208*p[1]/3 + 76*p[2] -112*p[3]/3 + 22*p[4]/3,
                  -80*p[0]/3 + 96*p[1] - 128*p[2] + 224*p[3]/3 - 16*p[4],
                  32*p[0]/3 - 128*p[1]/3 + 64*p[2] - 128*p[3]/3 + 32*p[4]/3
        };
    }

	template<typename Float, typename T>	
	constexpr T at(Float t, const std::array<T,samples>& p) const {
		auto c = coefficients(p);
		return (((c[4]*t + c[3])*t + c[2])*t + c[1])*t + c[0];	
	}
	
	template<typename Float, typename T>	
	constexpr T subrange(Float a, Float b, const std::array<T,samples>& p) const {
		auto c = coefficients(p);
		return ((((c[4]*b/5.0 + c[3]/4.0)*b + c[2]/3.0)*b + c[1]/2.0)*b + c[0])*b -
			   ((((c[4]*a/5.0 + c[3]/4.0)*a + c[2]/3.0)*a + c[1]/2.0)*a + c[0])*a;	
	}


/*    
	template<typename Float, typename T>	
	auto at(Float t, const std::array<T,samples>& p) const -> T {
		T c4 = 32*p[0]/3 - 128*p[1]/3 + 64*p[2] - 128*p[3]/3 + 32*p[4]/3;
		T c3 = -80*p[0]/3 + 96*p[1] - 128*p[2] + 224*p[3]/3 - 16*p[4];
		T c2 = 70*p[0]/3 -208*p[1]/3 + 76*p[2] -112*p[3]/3 + 22*p[4]/3;
		T c1 = -25*p[0]/3+ 16*p[1] - 12*p[2] +16*p[3]/3 - p[4];
		T c0 = p[0];
		return (((c4*t + c3)*t + c2)*t + c1)*t + c0;	
	}
	
	template<typename Float, typename T>	
	auto subrange(Float a, Float b, const std::array<T,samples>& p) const -> T {
		T c4 = 32*p[0]/3 - 128*p[1]/3 + 64*p[2] - 128*p[3]/3 + 32*p[4]/3;
		T c3 = -80*p[0]/3 + 96*p[1] - 128*p[2] + 224*p[3]/3 - 16*p[4];
		T c2 = 70*p[0]/3 -208*p[1]/3 + 76*p[2] -112*p[3]/3 + 22*p[4]/3;
		T c1 = -25*p[0]/3+ 16*p[1] - 12*p[2] +16*p[3]/3 - p[4];
		T c0 = p[0];
		return ((((c4*b/5.0 + c3/4.0)*b + c2/3.0)*b + c1/2.0)*b + c0)*b -
			   ((((c4*a/5.0 + c3/4.0)*a + c2/3.0)*a + c1/2.0)*a + c0)*a;	
	}
*/

} boole;

template<typename Q, std::size_t N>
struct Steps {
private:
	Q quadrature;
public:
	static constexpr std::size_t samples = (Q::samples - 1)*N + 1;

	Steps(const Q& quadrature) noexcept : quadrature(quadrature) {}
	
	template<typename T>
	auto operator()(const std::array<T,samples>& p) const -> T {
		std::array<T,Q::samples> local;
		std::copy(p.begin(),p.begin()+Q::samples,local.begin());
		T sol = quadrature(local)/N;
//		std::cerr<<std::endl<<0<<" -> "<<sol<<std::endl;
		for (std::size_t i = 1; i<N; ++i) {
			std::copy(p.begin()+i*(Q::samples - 1),
				      p.begin()+i*(Q::samples - 1) + Q::samples, 
					  local.begin());
			sol += quadrature(local)/N;
//			std::cerr<<i<<" -> "<<sol<<std::endl;
		}
		return sol;
	}
	
	template<typename Float, typename T>	
	auto at(Float t, const std::array<T,samples>& p) const -> T {
		std::size_t i = std::min(std::size_t(t*N),N-1);
		Float t_local = Float(t*N) - Float(i);
		std::array<T,Q::samples> local;
		std::copy(p.begin()+i*(Q::samples - 1),
			p.begin()+i*(Q::samples - 1) + Q::samples, 
			local.begin());
		return quadrature.at(t_local,local);
	}
	
	template<typename Float, typename T>
	auto subrange(Float a, Float b, const std::array<T,samples>& p) const -> T {
		std::size_t ia = std::min(std::size_t(a*N),N-1);
		std::size_t ib = std::min(std::size_t(b*N),N-1);
		Float a_local = Float(a*N) - Float(ia);
		Float b_local = Float(b*N) - Float(ib);
		std::array<T,Q::samples> local;
		std::copy(p.begin()+ia*(Q::samples - 1),
			p.begin()+ia*(Q::samples - 1) + Q::samples, 
			local.begin());
		if (ia == ib) return quadrature.subrange(a_local, b_local, local)/N;
		else {
			T sol = quadrature.subrange(a_local,Float(1), local)/N;
			for (std::size_t i = (ia+1); i<ib; ++i) {
				std::copy(p.begin()+i*(Q::samples - 1),
				      p.begin()+i*(Q::samples - 1) + Q::samples, 
					  local.begin());
				sol += quadrature(local)/N;
			}
			std::copy(p.begin()+ib*(Q::samples - 1),
				      p.begin()+ib*(Q::samples - 1) + Q::samples, 
					  local.begin());
			sol += quadrature.subrange(Float(0),b_local, local)/N;
			return sol;
		}
	}
};

template<std::size_t N, typename Q> 
auto steps(const Q& quadrature) {
	return Steps<Q,N>(quadrature);
}

}
