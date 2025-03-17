#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

class Config{ //Singleton class for configuration
public:
    static Config& getInstance();

private:
    bool m_SUBCYCLING = true; //advancing particles via subcycling
public:
    bool setSUBCYCLING(bool option);
    bool getSUBCYCLING() const;

private:
    bool m_MULTITHREADING = true; 
    unsigned int m_NUM_THREADS = 1;
public:
    bool setMULTITHREADING(bool option);
    bool getMULTITHREADING() const;
    unsigned int setNUM_THREADS(unsigned int n);
    unsigned int getNUM_THREADS() const;

private:
    bool m_MERGING = true; //merging particles Species::merge()
public:
    bool setMERGING(bool option);
    bool getMERGING() const;

private:
    bool m_INFLUENCE_ON_BACKGROUND = false;
public:
    bool setINFLUENCE_ON_BACKGROUND(bool option);
    bool getINFLUENCE_ON_BACKGROUND() const;

private:
    bool m_SPUTTERING = false;
public:
    bool setSPUTTERING(bool option);
    bool getSPUTTERING() const;


    Config(const Config& other) = delete;
    void operator=(const Config& other) = delete;

private:
    Config();
            
    // Output::modes OUTPUT = Output::modes::fields;
};

//////////////// functions ////////////////////////////////////

std::vector<size_t> splitIntoChunks(size_t size);




#endif // CONFIG_H