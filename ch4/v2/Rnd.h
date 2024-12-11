#ifndef RND_H
#define RND_H
#include <random>
#include "all.h"

class Rnd{
protected:
    std::mt19937 mt_gen;
    std::uniform_real_distribution<type_calc> rnd_dist;

public:
    Rnd(); // mt19937 has default seed;
    Rnd(unsigned seed); //our own seed
    type_calc operator()();
    type_calc operator()(type_calc min, type_calc max);
};

#endif