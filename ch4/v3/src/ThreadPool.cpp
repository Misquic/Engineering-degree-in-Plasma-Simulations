#include "ThreadPool.h"
#include "Config.h"

void ThreadWorker::operator()(){
    while(!thread_pool->m_stop || (thread_pool->m_stop && !thread_pool->m_tasks.empty())){
        std::unique_lock<std::mutex> lock(thread_pool->m_queue_mutex);

        if(thread_pool->m_num_tasks_pending == 0 && thread_pool->m_tasks.empty()){
            thread_pool->m_busy_cv.notify_all();
        }

        thread_pool->m_queue_cv.wait(lock, [this]{return this->thread_pool->m_stop || !this->thread_pool->m_tasks.empty();}); //thread waits until one of conditions is true
        if(!this->thread_pool->m_tasks.empty()){
#if __cpp_lib_move_only_function
            auto task = std::move(thread_pool->tasks.front());
#else
            auto task = thread_pool->m_tasks.front();
#endif
            thread_pool->m_tasks.pop();
            // thread_pool->m_num_tasks_pending++;
            lock.unlock();
            // using namespace std::chrono_literals;
            task();
            // std::this_thread::sleep_for(1ms);
            lock.lock();

            // thread_pool->m_busy_cv.notify_all();


            thread_pool->m_num_tasks_pending--;
            if(thread_pool->m_num_tasks_pending == 0 && thread_pool->m_tasks.empty()){
                thread_pool->m_busy_cv.notify_all();
            }
        }
    }
};

void ThreadPool::waitForCompletion(){
    while(true){
        std::unique_lock<std::mutex> lock_queue(m_queue_mutex);
        if(m_num_tasks_pending == 0 && m_tasks.empty()){
            return;
        }    
        lock_queue.unlock();
    }
    // std::unique_lock<std::mutex> lock_queue(m_queue_mutex);
    // if(m_num_tasks_pending == 0 && m_tasks.empty()){
    //     return;
    // }
    // lock_queue.unlock();
    // std::unique_lock<std::mutex> lock_busy(m_busy_mutex);
    
    // m_busy_cv.wait(lock_busy, [this]{std::unique_lock<std::mutex> lock_queue(m_queue_mutex);return m_num_tasks_pending == 0 && m_tasks.empty();});

    // m_busy_cv.wait(lock_busy, [this]{m_queue_cv.notify_all();std::unique_lock<std::mutex> lock(m_queue_mutex); std::cerr << "pending: " << m_num_tasks_pending << " queue: " << m_tasks.empty() << "\n" <<std::flush; return (m_num_tasks_pending == 0 && m_tasks.empty());});
    // std::cerr << "completed\n" << std::flush;
};


ThreadPool::ThreadPool(const unsigned int num_threads): m_num_tasks_pending(0), m_thread_workers(std::vector<std::thread>(num_threads)){
    for(size_t i = 0; i < num_threads; i++){
        m_thread_workers[i] = std::thread(ThreadWorker(this));
    }
};

void ThreadPool::shutDown(){
    {
        std::unique_lock<std::mutex> lock(m_queue_mutex);
        m_stop = true;
        m_queue_cv.notify_all();
    }
    for(std::thread& th: m_thread_workers){
        if(th.joinable()){
            th.join();
        }
    }
};

size_t ThreadPool::queueSize() const{
    std::unique_lock<std::mutex> lock(m_queue_mutex);
    return m_tasks.size();
};

size_t ThreadPool::numThreads()const{
    return m_thread_workers.size();
};

void ThreadPool::resize(unsigned int new_size){
    if(new_size > m_thread_workers.size()){
        m_thread_workers.reserve(new_size);
        for(unsigned int i = m_thread_workers.size(); i < new_size; i++){
            m_thread_workers.emplace_back(std::thread(ThreadWorker(this)));
        }
    }else if(new_size > 0 && new_size < m_thread_workers.size()){
        std::unique_lock<std::mutex> lock(m_queue_mutex);
        m_queue_cv.wait(lock, [this]{return m_num_tasks_pending == 0;});

        for(unsigned int i = m_thread_workers.size()-1; i >= new_size; i--){
            m_thread_workers[i].join();
            m_thread_workers.pop_back();
        }
    }
};


ThreadPool& SingletonThreadPool::getInstance(){
    static ThreadPool instance(Config::getInstance().getNUM_THREADS());

    unsigned int config_num_threads = Config::getInstance().getNUM_THREADS();
    if(instance.numThreads()!=config_num_threads){
        instance.resize(config_num_threads);
    }

    return instance;
};
