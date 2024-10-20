#include "Rnd.h"

/*constructors*/
Rnd::Rnd(): mt_gen{std::random_device()()}, rnd_dist{0.0, 1.0} { //[0, 1.o)
};
Rnd::Rnd(unsigned seed): mt_gen{seed}, rnd_dist{0.0, 1.0} { //[0, 1.o)
};

/*operators*/
type_calc Rnd::operator()(){
    return rnd_dist(mt_gen);
};
type_calc Rnd::operator()(type_calc min, type_calc max){
    return min + rnd_dist(mt_gen)*(max - min); //calculating that outside would be faster? less arguments to pass
};

Rnd rnd;