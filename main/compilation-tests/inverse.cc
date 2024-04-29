#include <math.h>
#include <iostream>
#include <array>
#include <vector>

float cuberoot(float x){ 
    float cube_root;
    if(x < 0)
    { 
        x = abs(x);
        cube_root = pow(x,1./3.)*(-1);
    }
    else{
        cube_root = pow(x,1./3.);
    }
    return cube_root;
}

std::array<float, 3> coefficients(const std::array<float,3>& p) {
        return std::array<float, 3>{
                  p[0],
                  -3*p[0]+4*p[1]-p[2],
                  2*p[0]-4*p[1]+2*p[2]
        };
    }

float cdf(float x, const std::array<float,3>& p)
{
    auto c = coefficients(p);
		return ((c[2]*x/3.0 + c[1]/2.0)*x + c[0])*x;
}

std::vector<float> inv_cdf(float x, const std::array<float,3>& points)
{
    std::vector<float> solutions;

    auto coeff = coefficients(points);
    float a = coeff[2]/3.;
    float b = coeff[1]/2.;
    float c = coeff[0];
    float d = -x;

    float p = c/a - pow(b,2.)/(3.*pow(a,2.));
    float q = 2*pow(b,3.)/(27.*pow(a,3.)) - b*c/(3.*pow(a,2.)) + d/a;

    float sqr = pow(q/2.,2.) + pow(p/3.,3.);
    if(sqr >= 0.){ 
        //One solution
        solutions.push_back(cuberoot(-q/2. - sqrt(sqr)) + cuberoot(-q/2. + sqrt(sqr))-b/(3*a));
    }
    else{ 
        //3 solutions
        for(int i=0; i<3; i++)
            solutions.push_back(2.*sqrt(-p/3.) * cos(1./3. * acos(3.*q/(2.*p) * sqrt(-3./p)) - 2.*M_PI * i/3.)-b/(3.*a));
    }
    return solutions;
}


/*Returns a sample from the inverted CDF normalized
    s -> Random uniform sample between 0-1
    a -> Minimus range value
    b -> Maximum range value
    p -> Polynomial sample points
*/
float sample_normalized(float s, float a, float b, const std::array<float,3>& p){ 
    float x = s*(cdf(b,p) - cdf(a,p)) + cdf(a,p);
    auto res = inv_cdf(x,p);
    bool solved = false; 
    float solution = 0;
    for(auto r : res){
        if(r > a && r < b){ 
            if(solved) std::cout<<"Warning, more than one solutions for range "<<a<<" - "<<b<<std::endl;
            solved = true;
            solution = r;
        }
    }
    if(!solved) std::cout<<"Warning, no solutions for range "<<a<<" - "<<b<<std::endl;
    return solution;
}

int main(int argc, char *argv[]){

    float x;
    std::array<float,3> p = {0.43,0.9,0.42};
    

    for(int i=0; i<10; i++)
        std::cout<<sample_normalized(i*0.1,0,0.125,p)<<std::endl;

    return 0;
}