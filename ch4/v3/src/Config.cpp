#include <thread>
#include "Config.h"
#include "all.h"

Config& Config::getInstance(){
    static Config instance;
    return instance;
};

bool Config::setSUBCYCLING(bool option){
    m_SUBCYCLING = option;
    return m_SUBCYCLING;
};
bool Config::getSUBCYCLING() const{
    return m_SUBCYCLING;
};

bool Config::setMULTITHREADING(bool option){
    m_MULTITHREADING = option;
    return m_MULTITHREADING;
};
bool Config::getMULTITHREADING() const{
    return m_MULTITHREADING;
};
unsigned int Config::setNUM_THREADS(unsigned int n){
    if(std::thread::hardware_concurrency() <= 1){ // if only one core avaiable
        dmsg("setNUM_THREADS not succesfull, only one core avaiable\n");
        m_NUM_THREADS = 1;
        m_MULTITHREADING = false;
    }else{
        if(n < 1 || n > std::thread::hardware_concurrency()){//allows for setting 1, for testing
            dmsg("setNUM_THREADS not succesfull, wrong value\n");
            return m_NUM_THREADS;
        }
        
        m_NUM_THREADS = n;
    }
    return m_NUM_THREADS;
};
unsigned int Config::getNUM_THREADS() const{
    return m_NUM_THREADS;
};

bool Config::setMERGING(bool option){
    m_MERGING = option;
    return m_MERGING;
};
bool Config::getMERGING() const{
    return m_MERGING;
};

bool Config::setINFLUENCE_ON_BACKGROUND(bool option){
    m_INFLUENCE_ON_BACKGROUND = option;
    return m_INFLUENCE_ON_BACKGROUND;
};
bool Config::getINFLUENCE_ON_BACKGROUND() const{
    return m_INFLUENCE_ON_BACKGROUND;
};

bool Config::setSPUTTERING(bool option){
    m_SPUTTERING = option;
    return m_SPUTTERING;
};
bool Config::getSPUTTERING() const{
    return m_SPUTTERING;
};

Config::Config(){
    if(std::thread::hardware_concurrency() <= 1){ // if only one core avaiable
        m_NUM_THREADS = 1;
        m_MULTITHREADING = false;
    }else{
        m_NUM_THREADS = std::thread::hardware_concurrency() - 1; //by default use all but one avaiable cores
    }
}; 

//////////////// functions ////////////////////////////////////

std::vector<size_t> splitIntoChunks(size_t size){// returns starting indexes of chunked vector/array starting index: indexes[thread_id], ending index (dont include): indexes[thread_id+1]
    //size - size of vector/array that needs to be splitted into chunks
    unsigned int num_threads = Config::getInstance().getNUM_THREADS();
    if(size < num_threads){
        std::vector<size_t> result(size+1);
        for(size_t i = 0; i < size; i++){
            result[i] = i;
        }
        result[size] = size;
        return result;
    }
    unsigned int more_than_base = size % num_threads;
    unsigned int base = size/num_threads;

    std::vector<size_t> indexes(num_threads+1);
    indexes[0] = 0;
    for(unsigned int i = 1; i <= num_threads; i++){
        indexes[i] = indexes[i-1] + base + (i<=more_than_base?1:0);
    }
    return indexes;
};


