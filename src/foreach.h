#pragma once
#include "multidimensional-range.h"
#include "log.h"
#include "progress-logger.h"
#include <thread>
#include <atomic>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <execution>
#include <cmath>

namespace viltrum {
inline struct _Sequential {} sequential;
inline struct _Parallel {} parallel;
inline struct _ParallelChunked {} parallel_chunked;

template<std::size_t DIM, typename F>
void for_each(const _Sequential& seq, const MultidimensionalRange<DIM>& range, const F& f) {
    for (auto pos : range) f(pos);
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const _Sequential& seq, const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    std::size_t total = range.size();   
    std::size_t r = 0;
    for (auto pos : range) {
        logger.log_progress(r++,total);
        f(pos);
    }
    logger.log_progress(total,total);
}

template<std::size_t DIM, typename F>
void for_each(const MultidimensionalRange<DIM>& range, const F& f) {
    for_each(sequential,range,f);
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    for_each(sequential,range,f,logger);
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const _Parallel& par, const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    std::size_t ntasks = range.max(0)-range.min(0); 
    std::vector<std::size_t> idxs(ntasks);
    std::atomic<std::size_t> done{0}; //Completed tasks
    std::iota(idxs.begin(), idxs.end(), 0);
    auto for_log = progress_logger(logger,ntasks,[&done] () { return done.load(std::memory_order_relaxed); });
    if constexpr (DIM > 1) {
        std::array<std::size_t,DIM-1> a, b;
        for (std::size_t s = 1;s<DIM;++s) {
            a[s-1] = range.min(s);
            b[s-1] = range.max(s);
        }
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) {
                for (auto pos : multidimensional_range(a,b)) f(d|pos);
                done.fetch_add(1,std::memory_order_relaxed);  
            });
    } else {
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) { 
                f(std::array<std::size_t,1>{d}); 
                done.fetch_add(1,std::memory_order_relaxed);  
            });
    }
    for_log.finish();
}

template<std::size_t DIM, typename F>
void for_each(const _Parallel& par, const MultidimensionalRange<DIM>& range, const F& f) {
    std::size_t ntasks = range.max(0)-range.min(0); 
    std::vector<std::size_t> idxs(ntasks);
    std::iota(idxs.begin(), idxs.end(), 0);
    if constexpr (DIM > 1) {
        std::array<std::size_t,DIM-1> a, b;
        for (std::size_t s = 1;s<DIM;++s) {
            a[s-1] = range.min(s);
            b[s-1] = range.max(s);
        }
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) {
                for (auto pos : multidimensional_range(a,b)) f(d|pos);
            });
    } else {
        std::for_each(std::execution::par_unseq,
            idxs.begin(),idxs.end(),
            [&] (std::size_t d) { 
                f(std::array<std::size_t,1>{d}); 
            });
    } 
}

template<std::size_t DIM, typename F, typename Logger>
void for_each(const _ParallelChunked& par, const MultidimensionalRange<DIM>& range, const F& f, Logger& logger) {
    // Determine chunk size per dimension based on hardware concurrency
    unsigned int nthreads = std::max(1u, std::thread::hardware_concurrency());
    unsigned int chunks_per_dim = std::max(1u, static_cast<unsigned int>(std::pow(nthreads, 1.0/DIM)));
    
    // Calculate chunk ranges for each dimension
    std::array<std::size_t, DIM> chunk_mins, chunk_maxs, chunk_sizes;
    std::size_t total_chunks = 1;
    for (std::size_t d = 0; d < DIM; ++d) {
        chunk_mins[d] = range.min(d);
        chunk_maxs[d] = range.max(d);
        std::size_t dim_size = chunk_maxs[d] - chunk_mins[d];
        chunk_sizes[d] = std::max(std::size_t(1), (dim_size + chunks_per_dim - 1) / chunks_per_dim);
        std::size_t num_chunks_in_dim = (dim_size + chunk_sizes[d] - 1) / chunk_sizes[d];
        total_chunks *= num_chunks_in_dim;
    }
    
    // Create chunk indices
    std::vector<std::size_t> chunk_idxs(total_chunks);
    std::iota(chunk_idxs.begin(), chunk_idxs.end(), 0);
    std::atomic<std::size_t> done{0}; //Completed chunks
    
    // Logger thread
    auto for_log = progress_logger(logger,total_chunks,[&done]() { return done.load(std::memory_order_relaxed); });
    
    // Process chunks in parallel
    std::for_each(std::execution::par_unseq,
        chunk_idxs.begin(), chunk_idxs.end(),
        [&](std::size_t chunk_idx) {
            // Convert linear chunk index to multidimensional chunk coordinates
            std::array<std::size_t, DIM> chunk_coord;
            std::size_t temp = chunk_idx;
            for (std::size_t d = 0; d < DIM; ++d) {
                std::size_t dim_size = chunk_maxs[d] - chunk_mins[d];
                std::size_t num_chunks_in_dim = (dim_size + chunk_sizes[d] - 1) / chunk_sizes[d];
                chunk_coord[d] = temp % num_chunks_in_dim;
                temp /= num_chunks_in_dim;
            }
            
            // Calculate actual range for this chunk
            std::array<std::size_t, DIM> chunk_start, chunk_end;
            for (std::size_t d = 0; d < DIM; ++d) {
                chunk_start[d] = chunk_mins[d] + chunk_coord[d] * chunk_sizes[d];
                chunk_end[d] = std::min(chunk_start[d] + chunk_sizes[d], chunk_maxs[d]);
            }
            
            // Process all positions in this chunk
            for (auto pos : multidimensional_range(chunk_start, chunk_end)) {
                f(pos);
            }
            done.fetch_add(1,std::memory_order_relaxed);
        });
    
    for_log.finish();
}

template<std::size_t DIM, typename F>
void for_each(const _ParallelChunked& par, const MultidimensionalRange<DIM>& range, const F& f) {
    // Determine chunk size per dimension based on hardware concurrency
    unsigned int nthreads = std::max(1u, std::thread::hardware_concurrency());
    unsigned int chunks_per_dim = std::max(1u, static_cast<unsigned int>(std::pow(nthreads, 1.0/DIM)));
    
    // Calculate chunk ranges for each dimension
    std::array<std::size_t, DIM> chunk_mins, chunk_maxs, chunk_sizes;
    std::size_t total_chunks = 1;
    for (std::size_t d = 0; d < DIM; ++d) {
        chunk_mins[d] = range.min(d);
        chunk_maxs[d] = range.max(d);
        std::size_t dim_size = chunk_maxs[d] - chunk_mins[d];
        chunk_sizes[d] = std::max(std::size_t(1), (dim_size + chunks_per_dim - 1) / chunks_per_dim);
        std::size_t num_chunks_in_dim = (dim_size + chunk_sizes[d] - 1) / chunk_sizes[d];
        total_chunks *= num_chunks_in_dim;
    }
    
    // Create chunk indices
    std::vector<std::size_t> chunk_idxs(total_chunks);
    std::iota(chunk_idxs.begin(), chunk_idxs.end(), 0);
    
    // Process chunks in parallel
    std::for_each(std::execution::par_unseq,
        chunk_idxs.begin(), chunk_idxs.end(),
        [&](std::size_t chunk_idx) {
            // Convert linear chunk index to multidimensional chunk coordinates
            std::array<std::size_t, DIM> chunk_coord;
            std::size_t temp = chunk_idx;
            for (std::size_t d = 0; d < DIM; ++d) {
                std::size_t dim_size = chunk_maxs[d] - chunk_mins[d];
                std::size_t num_chunks_in_dim = (dim_size + chunk_sizes[d] - 1) / chunk_sizes[d];
                chunk_coord[d] = temp % num_chunks_in_dim;
                temp /= num_chunks_in_dim;
            }
            
            // Calculate actual range for this chunk
            std::array<std::size_t, DIM> chunk_start, chunk_end;
            for (std::size_t d = 0; d < DIM; ++d) {
                chunk_start[d] = chunk_mins[d] + chunk_coord[d] * chunk_sizes[d];
                chunk_end[d] = std::min(chunk_start[d] + chunk_sizes[d], chunk_maxs[d]);
            }
            
            // Process all positions in this chunk
            for (auto pos : multidimensional_range(chunk_start, chunk_end)) {
                f(pos);
            }
        });
}

}