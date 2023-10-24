#include "../viltrum.h"
#include <iostream>

int main(int argc, char **argv) {
    for (auto x : viltrum::multidimensional_range(std::array<std::size_t,3>{1UL,2UL,3UL})) {
        std::cout<<"[ ";
        for (auto d : x) std::cout<<d<<" ";
        std::cout<<"]"<<std::endl;
    }
    std::cout<<std::endl;

    for (auto x : viltrum::multidimensional_range(std::array<std::size_t,2>{1UL,4UL}, std::array<std::size_t,2>{4UL,6UL})) {
        std::cout<<"[ ";
        for (auto d : x) std::cout<<d<<" ";
        std::cout<<"]"<<std::endl;
    }
}
