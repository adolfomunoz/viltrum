#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace viltrum {

/**
 * Accumulates, without any locking, the contributions that the processing of ONE region adds to the 
 * bins, and adds them to the shared bins (through a MutexedBins) in one go when flush() is called.
 *
 * Why: the contributions of a region (one for each sample of the residual, plus the integral of
 * the control variate in each bin) fall in the very few bins that the region overlaps, so with a direct
 * locked add there is one lock acquisition per sample, and all of them on a handful of locks. 
 * When the cost of each sample is small (for example with a cheap control variate), threads end up 
 * waiting for those locks and most of the time is spent in the mutexes, with a cost per add 
 * that is several times the cost of evaluating the sample. Accumulating locally reduces the 
 * number of locked adds from one per sample to one per bin overlapped by the region.
 * 
 * Implementation: the bins overlapped by the region (an axis aligned box of bins, which is what
 * pixels_in_region returns) are stored densely, so adding is just an index computation. 
 * Contributions that are outside of that box, or all of them if the box is too big to be worth
 * (more than max_dense bins), are directly added to the shared bins, so the result is always 
 * correct whatever the box is. 
 *
 * The sums are done in the order of the contributions of the region, and then the partial sums are 
 * added to the bins, so the result is equal to the direct add up to floating point rounding 
 * (the order of the additions to the bins was already random because of the parallelism).
 *
 * This class is not thread safe, it is meant to be created by each thread for each region.
 *
 * @tparam MutexedBins Type of the shared bins (MutexedBins or MutexedBinsAligned) with an add(pos,value) method.
 * @tparam Value The type of the values of the bins.
 */
template<typename MutexedBins, std::size_t DIMBINS, typename Value>
class RegionBinAccumulator {
    MutexedBins& shared;
    std::array<std::size_t,DIMBINS> start, extent, stride;
    std::vector<Value> values;
    std::vector<char> touched; //Only the touched bins are added to the shared ones when flushing

    //Returns false if the position is outside of the dense box
    bool index_of(const std::array<std::size_t,DIMBINS>& pos, std::size_t& index) const {
        index = 0;
        for (std::size_t d = 0; d<DIMBINS; ++d) {
            const std::size_t local = pos[d] - start[d]; //Wraps around (and is not less than extent) if pos[d] < start[d]
            if (local >= extent[d]) return false;
            index += local*stride[d];
        }
        return true;
    }

public:
    static constexpr std::size_t default_max_dense = 4096;

    /**
     * @param box The bins to accumulate densely: anything with min(d) and max(d) (max exclusive), 
     * like the MultidimensionalRange returned by pixels_in_region.
     */
    template<typename Box>
    RegionBinAccumulator(MutexedBins& shared, const Box& box, std::size_t max_dense = default_max_dense) : shared(shared) {
        std::size_t total = 1;
        bool dense = true;
        for (std::size_t d = 0; d<DIMBINS; ++d) {
            start[d] = box.min(d);
            extent[d] = (box.max(d) > box.min(d))?(box.max(d) - box.min(d)):0;
            stride[d] = total;
            if (extent[d] == 0 || total > max_dense/extent[d]) dense = false;
            else total *= extent[d];
        }
        if (dense) {
            values.assign(total,Value{});
            touched.assign(total,0);
        } else { //Everything goes directly to the shared bins
            extent.fill(0);
        }
    }

    RegionBinAccumulator(const RegionBinAccumulator&) = delete;
    RegionBinAccumulator& operator=(const RegionBinAccumulator&) = delete;

    template<typename V>
    void add(const std::array<std::size_t,DIMBINS>& pos, const V& v) {
        std::size_t index;
        if (index_of(pos,index)) {
            values[index] += v;
            touched[index] = 1;
        } else shared.add(pos,v);
    }

    //Adds all the accumulated contributions to the shared bins, and leaves the accumulator empty.
    void flush() {
        if (values.empty()) return;
        std::array<std::size_t,DIMBINS> pos = start;
        for (std::size_t index = 0; index<values.size(); ++index) {
            if (touched[index]) {
                shared.add(pos,values[index]);
                values[index] = Value{};
                touched[index] = 0;
            }
            //Next position of the box, with the first dimension varying the fastest (same as the strides)
            for (std::size_t d = 0; d<DIMBINS; ++d) {
                if (++pos[d] < (start[d] + extent[d])) break;
                pos[d] = start[d];
            }
        }
    }

    //Anything not flushed would be lost, so we do it automatically
    ~RegionBinAccumulator() { flush(); }
};

}
