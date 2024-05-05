#pragma once
#include <array>
#include <vector>
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

	template<typename Float, typename T>	
	constexpr T subrange(Float a, Float b, const std::array<T,samples>& p) const {
	    auto c = coefficients(p);
		return ((c[2]*b/3.0 + c[1]/2.0)*b + c[0])*b - ((c[2]*a/3.0 + c[1]/2.0)*a + c[0])*a;	
	}

	//My functions are below this line

	//I had to create this function, otherwise c++ returns -nans for a cube root of a negative number, even though 
	//it's a cube root
	constexpr static float cuberoot(float x) { 
		return (x<0)?(-1*pow(-x,1./3.)):pow(x,1./3.);
	}

	//Name cdf, inv_cdf and sample_normalized as you wish
	template<typename Float, typename T>
	constexpr T cdf(Float t, const std::array<T,samples>& p) const {
    	auto c = coefficients(p);
		return ((c[2]*t/3.0 + c[1]/2.0)*t + c[0])*t;
	}

	template<typename Float, typename T>
	std::vector<T> inv_cdf(Float x, const std::array<T,samples>& points) const {
		std::vector<T> solutions;

		auto coeff = coefficients(points);
		auto a = coeff[2]/3.;	//Term multiplying x³ in the cdf
		auto b = coeff[1]/2.;	//Term multiplying x² in the cdf
		auto c = coeff[0];		//Term multiplying x in the cdf
		auto d = -x;			//-x because we are solving the roots for cdf = x (cdf - x = 0) to compute the inverse

		if (a==0) { //Second degree equation	
			if (b==0) { //Degree one equation
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


	/*Returns a sample from the inverted CDF normalized
    s -> Random uniform sample between 0-1
    a -> Minimus range value
    b -> Maximum range value
    p -> Polynomial sample points
	*/
	template<typename Float, typename T, typename Norm = NormDefault>
	constexpr Float sample_normalized(Float s, Float a, Float b, const std::array<T,samples>& p,
			const Norm& norm = Norm()) const { 
		auto cdf_a = cdf(a,p);	//Just to save time (I don't know if compilers does this already)

		//Point where we are computing the inverse (s*(^F(b) - ^F(a)) + ^F(a))
		auto x = norm(s*(cdf(b,p) - cdf_a) + cdf_a);
		auto res = inv_cdf(x,p);

		//There are some warnings in two strange scenarious, hope we don't have to see them printed
		bool solved = false; 
		Float solution = s;
		for(auto r : res){
			if(r >= a && r <= b){ 
				if(solved) std::cout<<"Warning, more than one solutions for range "<<a<<" - "<<b<<std::endl;
				solved = true;
				solution = r;
			}
		}
		if(!solved) { 
			std::cout<<"Warning, no solutions for range "<<a<<" - "<<b;
			for (Float r : res) std::cout<<" "<<r;
			std::cout<<std::endl;
		}	
		return solution;
	}

	template<typename Float, typename T>	
	constexpr T pdf(Float t, const std::array<T,samples>& p, Float a, Float b) const {
		return at(t,p)/subrange(a,b,p);
	}


/*    
	template<typename Float, typename T>	
	auto at(Float t, const std::array<T,samples>& p) const -> T {
		T c2 = 2*p[0]-4*p[1]+2*p[2];
		T c1 = -3*p[0]+4*p[1]-p[2];
		T c0 = p[0];
		return (c2*t + c1)*t + c0;	
	}

	template<typename Float, typename T>	
	auto subrange(Float a, Float b, const std::array<T,samples>& p) const -> T {
		T c2 = 2*p[0]-4*p[1]+2*p[2];
		T c1 = -3*p[0]+4*p[1]-p[2];
		T c0 = p[0];
		return ((c2*b/3.0 + c1/2.0)*b + c0)*b - ((c2*a/3.0 + c1/2.0)*a + c0)*a;	
	}
*/

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
