#include "../../viltrum.h"

//We cannot use this function below in viltrum because we cannot know the data type of Seq.
/*
template<typename Seq>
float taylor_exp(const Seq& seq) {
    float sum = 0.0f; float term = 1.0f; unsigned int n = 0;    
    auto it = seq.begin(); float x = *it; ++it;
    while ((*it) < term) {
        sum += 1.0f; // sum = sum + term / term (because we use it as a probablility as well)
        ++it; ++n;
        term *= x/float(n); //The next term is the previous one times x/n, so we can update it iteratively. We also use it as a probability, so we don't need to divide by n!
    } 
    return sum;
}
*/

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
    unsigned long samples = 1024;
    for (int i = 0; i<argc-1; ++i) {
        if (std::string(argv[i])=="-samples") samples = atol(argv[++i]);
    }
    
    auto method = viltrum::monte_carlo(samples);
    auto range = viltrum::range_primary_infinite<float>();

    std::cout<<"exp integrated between 0 and 1 approx = 1.7182818"<<std::endl;
    std::cout<<viltrum::integrate(method,TaylorExp(),range)<<std::endl;
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

