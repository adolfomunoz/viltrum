/**
 * Developed within the framework of Diego Viñuales' final bachelor's degree project at the University of Zaragoza,
 * under the supervision of Dr. Adolfo Muñoz
 * 
 * Then it was modified to be more general.
 */

#pragma once
#include "nested.h"
#include "error-heuristic.h"
#include <thread>
#include <atomic>
#include <vector>
#include <numeric>
#include <execution>
#include <algorithm>
#include "../priorityqueues/lockedpriorityqueue.hpp"
#include "../priorityqueues/multiqueue.hpp"
#include "../priorityqueues/skiplist.hpp"




namespace viltrum {

  
struct Comparator {
    template<typename T>
    bool operator()(const T& a, const T& b) const {
        auto errA = std::get<0>(a.extra());
        auto errB = std::get<0>(b.extra());
        
        bool nanA = std::isnan(errA);
        bool nanB = std::isnan(errB);
        
        if (nanA && nanB) return false; 
        
        if (nanA) return true;  
        if (nanB) return false;
        
        return errA < errB;
    }
};

namespace heap {
    struct multiqueue {
        template<typename T, typename Comparator>
        static auto make(const Comparator& comp) {
           unsigned int num_threads = std::thread::hardware_concurrency();
           if (num_threads == 0) num_threads = 2;   
           return Multiqueue<T, Comparator>(4 * num_threads, 2, comp);
        }
    };
    struct lockedpriorityqueue {
        template<typename T, typename Comparator>
        static auto make(const Comparator& comp) {
            return LockedPriorityQueue<T, Comparator>(comp);
        }
    };
    // I am not adding the skip list because our experiments show that
    // there is a small chance of deadlock, and it actually
    // performs worse than the multiqueue in our experiments.
}


template<typename Rule, typename ErrorHeuristic, typename Structure = heap::multiqueue,
         typename = std::enable_if_t<is_nested<Rule>::value>>
class RegionsGeneratorAdaptiveParallel {
    Rule rule;
    ErrorHeuristic error_heuristic;
    std::size_t subdivisions;

public:
    RegionsGeneratorAdaptiveParallel(const Rule& r, const ErrorHeuristic& er, std::size_t subdivisions)
        : rule(r), error_heuristic(er), subdivisions(subdivisions) {}

    template<std::size_t DIMBINS, typename F, typename Float, std::size_t DIM, typename Logger>
    auto generate(const std::array<std::size_t, DIMBINS>& bin_resolution,
                  const F& f, const Range<Float, DIM>& range, Logger& logger) const {

        auto r_init = region(f, rule, range.min(), range.max());
        auto errdim_init = error_heuristic(r_init);
        using ERegion = ExtendedRegion<decltype(r_init), decltype(errdim_init)>;

        unsigned int num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) num_threads = 2;

        auto mq = Structure::template make<ERegion>(Comparator());
        mq.push(ERegion(r_init, errdim_init));

        std::atomic<std::size_t> completed{0};

        logger.log_progress(std::size_t(0), subdivisions);

        std::vector<unsigned int> workers(num_threads);
        std::iota(workers.begin(), workers.end(), 0u);

        // par (not par_unseq) because the workers take locks.
        std::for_each(std::execution::par, workers.begin(), workers.end(), [&](unsigned int worker) {
            while (true) {
                std::size_t current = completed.load(std::memory_order_relaxed);
                if (current >= subdivisions) break;

                std::optional<ERegion> opt_r = mq.try_pop();
                if (!opt_r) {
                    std::this_thread::yield();
                    continue;
                }

                std::size_t my_slot = completed.fetch_add(1, std::memory_order_relaxed);
                if (my_slot >= subdivisions) {
                    mq.push(*opt_r);
                    break;
                }

                auto subregions = opt_r->split(f, std::get<1>(opt_r->extra()));
                for (auto& sr : subregions)
                    mq.push(ERegion(sr, error_heuristic(sr)));

                // Only worker 0 logs, so the logger is never called concurrently.
                if ((worker == 0) && (completed.load(std::memory_order_relaxed) < subdivisions))
                    logger.log_progress(std::min(completed.load(std::memory_order_relaxed), subdivisions), subdivisions);
            }
        });

        std::vector<ERegion> final_heap = mq.drain();

        logger.log_progress(subdivisions, subdivisions);

        return final_heap;
    }
};

template<typename Rule, typename ErrorHeuristic, typename Structure>
auto regions_generator_adaptive_parallel(const Rule& r, const ErrorHeuristic& er, const Structure& structure,
                                           std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveParallel<Rule, ErrorHeuristic,Structure>(r, er, subdivisions);
}

template<typename Rule, typename ErrorHeuristic>
auto regions_generator_adaptive_parallel(const Rule& r, const ErrorHeuristic& er,
                                           std::size_t subdivisions) {
    return RegionsGeneratorAdaptiveParallel<Rule, ErrorHeuristic,heap::multiqueue>(r, er, subdivisions);
}

}