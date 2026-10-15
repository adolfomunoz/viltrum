/**
 * Developed within the framework of Diego Viñuales' final bachelor's degree project at the University of Zaragoza,
 * under the supervision of Dr. Adolfo Muñoz
 * 
 * From repository:
 * Relaxed Priority Queues - https://github.com/DiegoV8/relaxedPriorityQueues_874280.git
 */

#pragma once

#include <functional>
#include <mutex>
#include <optional>
#include <queue>
#include <vector>

// std::priority_queue protected by a single mutex. Compare orders the elements
// like in std::priority_queue (the "greatest" element according to it is popped first).
template <typename T, typename Compare = std::less<T>>
class LockedPriorityQueue {
    // Exposes the protected underlying container of std::priority_queue.
    struct Heap : std::priority_queue<T, std::vector<T>, Compare> {
        explicit Heap(const Compare& comp) : std::priority_queue<T, std::vector<T>, Compare>(comp) {}
        std::vector<T>& container() { return this->c; }
    };
    Heap pq;
    mutable std::mutex mtx;
public:
    using value_type = T;

    explicit LockedPriorityQueue(const Compare& comp = Compare()) : pq(comp) {}

    void push(T val) {
        std::lock_guard<std::mutex> lock(mtx);
        pq.push(val);
    }

    std::optional<T> try_pop() {
        std::lock_guard<std::mutex> lock(mtx);
        if (pq.empty()) return std::nullopt;
        T val = pq.top();
        pq.pop();
        return val;
    }

    bool empty() const {
        std::lock_guard<std::mutex> lock(mtx);
        return pq.empty();
    }

    std::size_t size() const {
        std::lock_guard<std::mutex> lock(mtx);
        return pq.size();
    }

    // Moves out all elements (unsorted), leaving the queue empty. Not thread-safe.
    std::vector<T> drain() {
        return std::move(pq.container());
    }
};