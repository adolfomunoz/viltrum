#include "../../src/newton-cotes/region.h"
#include <random>
#include <string>
#include <iomanip>

int main(int argc, char** argv) {
    float p0 = 1, p1 = 0.2, p2 = 1;
    std::size_t bins = 6, samples = 100;
    float xmin = 0.0, xmax = 1.0;

    for (int i = 0;i<(argc-1);++i) {
        if (std::string(argv[i])=="-bins") bins = std::atoi(argv[++i]);
        else if (std::string(argv[i])=="-samples") samples = std::atoi(argv[++i]);
        else if (std::string(argv[i])=="-p0") p0 = std::atof(argv[++i]);
        else if (std::string(argv[i])=="-p1") p1 = std::atof(argv[++i]);
        else if (std::string(argv[i])=="-p2") p2 = std::atof(argv[++i]);
        else if (std::string(argv[i])=="-xmin") xmin = std::atof(argv[++i]);
        else if (std::string(argv[i])=="-xmax") xmax = std::atof(argv[++i]);
    }

    auto f = [p0,p1,p2] (std::array<float,1> x) {
        return p0+(-3*p0+4*p1-p2)*x[0] + (2*p0-4*p1+2*p2)*x[0]*x[0];
    };

    auto r = viltrum::region(f,viltrum::simpson, viltrum::range(0.0f,1.0f));

    std::vector<float> probs(bins,0), hist(bins,0);

    float dx = (xmax-xmin)/bins; 
    for (std::size_t i = 0; i<bins; ++i) {
        std::array<float,1> x;
        x[0] = (float(i)+0.5)*dx + xmin;
        probs[i] = r.pdf_subrange(x,viltrum::range(xmin,xmax));
    }

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<float> dist(0,1);
    for (std::size_t i = 0; i<samples; ++i) {
        std::array<float,1> smp;
        smp[0]=dist(rng);
        auto x = r.sample_subrange(smp,viltrum::range(xmin,xmax));
        std::size_t bin = bins*(x[0]-xmin)/(xmax-xmin);
        hist[bin] += float(bins)/float(samples);
    }

    float tot = 0;
    for (float x : probs) { tot+=x; std::cout<<std::setw(6)<<std::setprecision(3)<<x<<" "; }
    std::cout<<" [ "<<std::setw(6)<<std::setprecision(3)<<tot<< " ] "<<std::endl;
    tot = 0;
    for (float x : hist) { tot+=x; std::cout<<std::setw(6)<<std::setprecision(3)<<x<<" "; }
    std::cout<<" [ "<<std::setw(6)<<std::setprecision(3)<<tot<< " ] "<<std::endl;
}