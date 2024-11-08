#ifndef SOURCE_H
#define SOURCE_H
#include "Species.h"
#include "World.h"
#include "Vec3.h"
#include "all.h"

/*ideas:
    make it choose a face of inlet, but overload or sth so its not computing consuming*/

class ColdBeamSource{ //Z- FACE
protected:
    Species   &sp;      //reference to injected species
    World     &world;   //reference to world
    type_calc v_drift;  //mean drift velocity
    type_calc den;      // injection density
    std::string inlet_Face;

    //precalc:
    type_calc3 dx;
    type_calc3 x0;
    type_calc Lx; //length x
    type_calc Ly; //length y
    type_calc A; //surface area
    type_calc num_micro; //number of microparticles



public:
    /*constructor*/
    ColdBeamSource(Species& species, World& world, type_calc v_drift, type_calc den, std::string inlet_Face = "-z") noexcept;

    /*methods*/
    void sample();

};

#endif