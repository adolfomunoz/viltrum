#pragma once

#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>

namespace viltrum {

class LoggerNull {
public:
    template<typename Number>
    void log_progress(const Number& number, const Number& last = Number(1)) {}
    template<typename Data>
    void log(const Data& data) {}
};

class LoggerProgress {
    std::string name;
    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now(); 
public:
    LoggerProgress(const std::string& n) : name(n) {}
    template<typename Number>
    void log_progress(const Number& number, const Number& last = Number(1)) {
        static Number prev_number(0);
        if (number<=Number(0)) {
            start = std::chrono::steady_clock::now();
            std::cerr<<name<<" -                           \r";
        } else if (std::abs(double(number)-double(prev_number))>(double(last)*0.01)) {
            auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::steady_clock::now() - start);
            std::cerr<<name<<" - \t"<<std::fixed<<std::setprecision(2)<<std::setw(6)<<(100.0f*float(number)/float(last))<<"%\t("<<std::setprecision(3)<<std::setw(6)<<elapsed.count()<<" seconds)\r";
            prev_number=number;
        }  
    }
    template<typename Data>
    void log(const Data& data) {}
};

};

