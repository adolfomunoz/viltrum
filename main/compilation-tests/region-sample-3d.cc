#include "../../src/newton-cotes/region.h"
#include "../../src/multidimensional-range.h"
#include <random>
#include <string>
#include <iomanip>
#include <cmath>

int main(int argc, char** argv) {
    std::size_t bins = 3, samples = 100;
    std::array<float,3> xmin{0,0,0}, xmax{1,1,1};

    for (int i = 0;i<(argc-1);++i) {
        if (std::string(argv[i])=="-bins") bins = std::atoi(argv[++i]);
        else if (std::string(argv[i])=="-samples") samples = std::atoi(argv[++i]);
    }
    for (int i = 0;i<(argc-3);++i) {
        if (std::string(argv[i])=="-xmin") { 
            xmin[0] = std::atof(argv[++i]);
            xmin[1] = std::atof(argv[++i]);
            xmin[2] = std::atof(argv[++i]);
        }
        else if (std::string(argv[i])=="-xmax") {
            xmax[0] = std::atof(argv[++i]);
            xmax[1] = std::atof(argv[++i]);
            xmax[2] = std::atof(argv[++i]);
        }
    }

    auto f = [] (std::array<float,3> x) {
        return std::sin(0.5*x[0]*M_PI) + std::cos(0.5*x[1]*M_PI) + x[2];
    };

    auto r = viltrum::region(f,viltrum::simpson, viltrum::range_primary<3>());

    std::vector<float> probs(bins*bins*bins,0), hist(bins*bins*bins,0);

    std::array<float,3> dx;
    for (std::size_t c = 0; c<3; ++c) dx[c] = (xmax[c]-xmin[c])/bins;
    for (auto p : viltrum::multidimensional_range(std::array<size_t,3>{bins,bins,bins})) {
        std::array<float,3> x;
        for (std::size_t c = 0; c<3; ++c) x[c] = (float(p[c])+0.5)*dx[c] + xmin[c];
        probs[p[0]+bins*p[1]+bins*bins*p[2]] = r.pdf_subrange(x,viltrum::range(xmin,xmax));
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<float> dist(0,1);
    for (std::size_t i = 0; i<samples; ++i) {
        std::array<float,3> smp;
        for (std::size_t c = 0; c<3; ++c) smp[c]=dist(rng);
        auto x = r.sample_subrange(smp,viltrum::range(xmin,xmax));
        std::array<std::size_t,3> p;
        for (std::size_t c = 0; c<3; ++c) p[c] = std::min(bins-1,std::size_t(bins*(x[c]-xmin[c])/(xmax[c]-xmin[c])));
        hist[p[0]+bins*p[1]+bins*bins*p[2]] += float(bins*bins*bins)/float(samples);
    }

    float tot = 0; std::size_t i = 0;

    for (float x : probs) { 
        tot+=x; 
        std::cout<<std::setw(6)<<std::setprecision(3)<<x<<" ";
        i++;
        if ((i%(bins*bins)) == 0) std::cout<<std::endl;
        else if ((i%bins) == 0) std::cout<<"| ";
    }
    std::cout<<" [ "<<std::setw(6)<<std::setprecision(3)<<tot<< " ] "<<std::endl;
    tot = 0; i = 0;
    for (float x : hist) { 
        tot+=x; 
        std::cout<<std::setw(6)<<std::setprecision(3)<<x<<" "; 
        i++;
        if ((i%(bins*bins)) == 0) std::cout<<std::endl;
        else if ((i%bins) == 0) std::cout<<"| ";
    }
    std::cout<<" [ "<<std::setw(6)<<std::setprecision(3)<<tot<< " ] "<<std::endl;
}