#pragma once
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>

namespace viltrum {

/**
 * Reports the progress of a parallel computation to a logger from a helper thread, while the 
 * worker threads only have to update a counter.
 *
 * Usage:
 *     std::atomic<std::size_t> progress{0};
 *     auto progress_logger = viltrum::progress_logger(logger, total, [&] { return progress.load(); });
 *     ... workers do ++progress ...
 *     progress_logger.finish();
 *
 * - The constructor logs the initial progress synchronously, and then starts a thread that 
 *   logs the progress every "period" (250 ms by default).
 * - The helper thread sleeps on a condition variable, so it uses no CPU while waiting and, 
 *   more importantly, finish() wakes it up immediately. With a plain sleep, finishing would
 *   wait for the rest of the period, which adds up to "period" of latency to every parallel 
 *   computation (and hides the real time of the ones that are shorter than that).
 * - finish() stops and joins the helper thread and logs the final progress (total, total).
 * - If finish() is never called (for example because an exception is thrown by the workers) the
 *   destructor still stops and joins the thread, but does not report the computation as completed.
 *
 * Progress is a callable without arguments that returns the amount of completed work 
 * (in the same units as total). It is called from the helper thread, so it should be thread safe 
 * (reading an atomic counter, for example). 
 * The logger must outlive the ProgressLogger, and it is only used from the helper thread 
 * while it is running (and from the thread that calls the constructor and finish, when the 
 * helper thread is not running or has not been started yet).
 */
template<typename Logger, typename Progress>
class ProgressLogger {
    Logger& logger;
    Progress progress;
    std::size_t total;
    std::chrono::milliseconds period;

    std::mutex mutex;
    std::condition_variable wake_up;
    bool stop_requested = false;
    std::thread helper;

    void run() {
        std::unique_lock<std::mutex> lock(mutex);
        // wait_for returns true as soon as the stop is requested, and false when the period elapsed.
        while (!wake_up.wait_for(lock,period,[this] { return stop_requested; })) {
            // Not holding the lock while logging, so stop() is never blocked by slow logging.
            lock.unlock();
            logger.log_progress(std::size_t(progress()),total);
            lock.lock();
        }
    }

    void stop() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            stop_requested = true;
        }
        wake_up.notify_one();
        if (helper.joinable()) helper.join();
    }

public:
    ProgressLogger(Logger& logger, std::size_t total, Progress progress, 
            std::chrono::milliseconds period = std::chrono::milliseconds(250)) 
        : logger(logger), progress(std::move(progress)), total(total), period(period) {
        this->logger.log_progress(std::size_t(this->progress()),this->total);
        helper = std::thread([this] { run(); });
    }

    //The helper thread refers to this object, so it cannot be copied or moved.
    ProgressLogger(const ProgressLogger&) = delete;
    ProgressLogger& operator=(const ProgressLogger&) = delete;

    ~ProgressLogger() { stop(); }

    //Stops the helper thread and logs that all the work is completed.
    void finish() {
        if (!helper.joinable()) return; //Already finished
        stop();
        logger.log_progress(total,total);
    }
};

template<typename Logger, typename Progress>
auto progress_logger(Logger& logger, std::size_t total, Progress&& progress,
        std::chrono::milliseconds period = std::chrono::milliseconds(250)) {
    return ProgressLogger<Logger,std::decay_t<Progress>>(logger,total,std::forward<Progress>(progress),period);
}

}
