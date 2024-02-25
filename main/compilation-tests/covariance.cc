#include "../../src/control-variates/weight-strategy.h"
#include <random>
#include <list>
#include <string>
#include <iostream>

//This tests the accuracy of the online calculation of variance and covariance of the
//cv_optimize_weight
using namespace std;
int main(int argc, const char** argv) {
    float randomness = 0.1;
    int nsamples = 20;
    for (int i = 0; i<argc-1; ++i) {
        if (string(argv[i]) == "-randomness") randomness = atof(argv[++i]);
        else if (string(argv[i]) == "-nsamples") nsamples = atoi(argv[++i]);
    }

    std::list<float> x0, x1;
    auto acc = viltrum::cv_optimize_weight().accumulator<float>();  

    uniform_real_distribution<float> rnd(-randomness,randomness);
    uniform_real_distribution<float> x(0,2*M_PI);
    random_device rd;
    mt19937 rng(rd());
    for (int i = 0; i<nsamples; ++i) {
        float v0 = std::sin(x(rng));
        float v1 = v0 + rnd(rng);
        x0.push_back(v0); x1.push_back(v1);
        acc.push(v1,v0);
    }

    float mean0=0;
    for (float v: x0) mean0+=v;
    mean0/=nsamples;
    float mean1=0;
    for (float v: x1) mean1+=v;
    mean1/=nsamples;

    float var0=0;
    for (float v: x0) var0+=(v-mean0)*(v-mean0); 
    var0/=nsamples;

    float co = 0;
    list<float>::const_iterator i0,i1;
    for (i0=x0.begin(),i1=x1.begin();(i0!=x0.end())&&(i1!=x1.end());++i0,++i1)
        co+=(*i0 - mean0)*(*i1 - mean1);
    co/=nsamples;

    cout<<"Variance   = "<<var0<<" vs "<<acc.variance()<<endl;
    cout<<"Covariance = "<< co <<" vs "<<acc.covariance()<<endl;  
} 