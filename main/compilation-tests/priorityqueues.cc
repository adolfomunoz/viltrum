/**
 * Developed within the framework of Diego Viñuales' final bachelor's degree 
 * project at the University of Zaragoza, supervised by Dr. Adolfo Muñoz
 * 
 * This file contains a benchmark for the priority queue implementations 
 * used in the project.
 * It measures the performance of concurrent insertion and extraction 
 * of elements.
 * 
 * This has been adapted / vibecoded from the original benchmark in the 
 * repository:
 * https://github.com/DiegoV8/testsEstructuraDatos_874280.git
 */


#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <functional>
#include <iostream>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <string>
#include <vector>
#include "../../src/priorityqueues/multiqueue.hpp"
#include "../../src/priorityqueues/skiplist.hpp"
#include "../../src/priorityqueues/lockedpriorityqueue.hpp"



// Runs f(task_index) for 'threads' parallel tasks (scheduled by TBB) and waits for all of them.
// std::execution::par (not par_unseq) because the tasks take locks.
template <typename F>
void run_parallel(int threads, F f) {
    std::vector<int> ids(threads);
    std::iota(ids.begin(), ids.end(), 0);
    std::for_each(std::execution::par, ids.begin(), ids.end(), f);
}

// Milliseconds taken by f().
template <typename F>
double time_ms(F f) {
    auto start = std::chrono::steady_clock::now();
    f();
    std::chrono::duration<double, std::milli> elapsed = std::chrono::steady_clock::now() - start;
    return elapsed.count();
}

template <typename Queue>
void run_benchmark(const std::string& name, Queue& q, int n, int threads, int show_first) {
    std::vector<int> data(n);
    std::iota(data.begin(), data.end(), 1);
    std::shuffle(data.begin(), data.end(), std::mt19937(std::random_device{}()));

    auto push_all = [&](int t) {
        for (int i = t; i < n; i += threads) q.push(data[i]);
    };

    // Pops exactly n elements; a failed try_pop is retried, so empty() is never needed.
    auto pop_all = [&] {
        std::atomic<int> popped{0};
        run_parallel(threads, [&](int) {
            while (popped.load() < n)
                if (q.try_pop()) ++popped;
        });
    };

    // Removes whatever is left (untimed cleanup between tests).
    auto clear = [&] {
        while (q.size() > 0) q.try_pop();
    };

    // Speed test: concurrent insertion followed by concurrent extraction.
    double insert_time = time_ms([&] { run_parallel(threads, push_all); });
    double extract_time = time_ms(pop_all);

    // Drain test: get all elements at once in a vector.
    run_parallel(threads, push_all);
    std::vector<int> drained;
    double drain_time = time_ms([&] { drained = q.drain(); });
    clear();

    // Mixed test: the queue holds n elements and each task alternates a push and a pop.
    run_parallel(threads, push_all);
    double mixed_time = time_ms([&] {
        run_parallel(threads, [&](int t) {
            for (int i = t; i < n; i += threads) {
                q.push(data[i]);
                q.try_pop();
            }
        });
    });
    clear();

    // Quality test: concurrent insertion, sequential extraction.
    run_parallel(threads, push_all);
    std::vector<int> extracted;
    while ((int)extracted.size() < n)
        if (auto v = q.try_pop()) extracted.push_back(*v);

    // The ideal extraction order is n, n-1, ..., 1.
    int inversions = 0;
    double error_sum = 0.0;
    for (size_t i = 0; i < extracted.size(); ++i) {
        if (i + 1 < extracted.size() && extracted[i] < extracted[i + 1]) ++inversions;
        error_sum += std::abs(extracted[i] - (n - (int)i));
    }
    double mean_error = extracted.empty() ? 0.0 : error_sum / extracted.size();
    double variance = 0.0;
    for (size_t i = 0; i < extracted.size(); ++i) {
        double e = std::abs(extracted[i] - (n - (int)i));
        variance += (e - mean_error) * (e - mean_error);
    }
    double std_dev = extracted.empty() ? 0.0 : std::sqrt(variance / extracted.size());

    std::cout << "=== " << name << " (" << threads << " threads, " << n << " elements) ===\n"
              << "Insertion time:       " << insert_time << " ms\n"
              << "Extraction time:      " << extract_time << " ms\n"
              << "Drain time:           " << drain_time << " ms (" << drained.size() << " elements)\n"
              << "Mixed push/pop time:  " << mixed_time << " ms (" << n << " push+pop pairs)\n"
              << "Extracted elements:   " << extracted.size() << " of " << n << "\n"
              << "Order inversions:     " << inversions << "\n"
              << "Mean position error:  " << mean_error << "\n"
              << "Standard deviation:   " << std_dev << "\n"
              << "First " << show_first << " extracted: ";
    for (size_t i = 0; i < extracted.size() && (int)i < show_first; ++i)
        std::cout << (i ? ", " : "") << extracted[i];
    std::cout << "\n\n";
}

int main(int argc, char** argv) {
    int n = argc > 1 ? std::atoi(argv[1]) : 100000;
    int threads = argc > 2 ? std::atoi(argv[2]) : 4;
    int show_first = argc > 3 ? std::atoi(argv[3]) : 10;
    if (n <= 0 || threads <= 0 || show_first < 0) {
        std::cerr << "Usage: " << argv[0] << " [elements] [threads] [shown_elements]\n";
        return 1;
    }

    LockedPriorityQueue<int> locked;
    run_benchmark("Mutex-protected std::priority_queue", locked, n, threads, show_first);

    Multiqueue<int, std::less<int>> multiqueue(2 * threads, 2, std::less<int>());
    run_benchmark("Multiqueue", multiqueue, n, threads, show_first);

    Skiplist<int, std::greater<int>> skiplist(std::greater<int>(), 2 * threads);
    run_benchmark("Skiplist", skiplist, n, threads, show_first);

    return 0;
}
