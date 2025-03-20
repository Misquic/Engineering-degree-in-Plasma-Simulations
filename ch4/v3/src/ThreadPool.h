#ifndef THREADPOOL_H
#define THREADPOOL_H
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <condition_variable>
#include <iostream>

class ThreadPool;

class ThreadWorker{
public:
    ThreadWorker(ThreadPool* pool): thread_pool(pool){};
    void operator()();
private:
        ThreadPool* thread_pool = nullptr;
};

class ThreadPool{
public:
    ThreadPool(const unsigned int num_threads = std::thread::hardware_concurrency() - 1);
    ~ThreadPool(){shutDown();};
    
    void shutDown();
    size_t queueSize() const;
    size_t numThreads() const;
    void resize(unsigned int new_size);
    void waitForCompletion();

    template <typename F, typename... Args>
    auto AddTask(F&& f, Args&&... args) -> std::future<decltype(f(args...))>;


    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&)      = delete;

    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&&)      = delete;

private:
    friend class ThreadWorker;
    // unsigned int m_busy_threads;
    std::atomic<unsigned int> m_num_tasks_pending;
    mutable std::mutex m_queue_mutex;
    mutable std::mutex m_busy_mutex;
    std::vector<std::thread> m_thread_workers;
#if __cpp_lib_move_only_function
    std::queue<std::move_only_function<void()>> m_tasks;
#else
    std::queue<std::function<void()>> m_tasks; //we want function that returns nothing, and takes no arguments
#endif
    std::condition_variable m_queue_cv;
    std::condition_variable m_busy_cv;
    bool m_stop = false;
};

template <typename F, typename... Args>
auto ThreadPool::AddTask(F&& f, Args&&... args) -> std::future<decltype(f(args...))>{
    m_num_tasks_pending++;
#if __cpp_lib_move_only_function
    std::packaged_task<decltype(f(args...))()> task(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    auto future       = task.get_future();
    auto wrapper_func = [task = move(task)]() mutable { std::move(task)(); };
    {
        std::lock_guard<std::mutex> lock(m_queue_mutex);
        m_tasks.push(std::move(wrapper_func));
        // Wake up one thread if its waiting
        m_queue_cv.notify_one();
    }

    // Return future from promise
    return future;
#else

    auto task_ptr = std::make_shared<std::packaged_task<decltype(f(args...))()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    auto wrapper_func = [task_ptr]() { (*task_ptr)(); };
    {
        std::lock_guard<std::mutex> lock(m_queue_mutex);
        m_tasks.push(wrapper_func);
        // Wake up one thread if its waiting
        m_queue_cv.notify_one();
    }

    // Return future from promise
    return task_ptr->get_future();
#endif
}

// template <typename F, typename... Args>
// void ThreadPool::AddTaskVoid(F&& f, Args&&... args){

//     //std::bind seems to be enaugh to make callable object that we want to call
// #if __cpp_lib_move_only_function

// #else    
//     std::bind(std::forward<F>(f), std::forward<Args>(args)...)

//     auto task_ptr = std::make_shared<std::packaged_task<decltype(f(args...))()>> (std::bind(std::forward<F>(f), std::forward<Args>(args)...));
//     auto wrapper_func = [task_ptr]() {(*task_ptr)();};
//     {
//         std::lock_guard<std::mutex> lock(m_queue_mutex);
//         m_tasks.push(wrapper_func);
//         m_queue_cv.notify_one();
//     }
//     return task_ptr->get_future();
// #endif
// };

class SingletonThreadPool{
    public:
        static ThreadPool& getInstance();
    
        SingletonThreadPool(const SingletonThreadPool& other) = delete;
        SingletonThreadPool(SingletonThreadPool&& other) = delete;
        void operator=(const SingletonThreadPool& other) = delete;
        void operator=(SingletonThreadPool&& other) = delete;
    
    private:
        static unsigned int m_num_threads;
        SingletonThreadPool() = default;
    };


#endif //THREADPOOL_H