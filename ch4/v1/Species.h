#ifndef SPECIES_H
#define SPECIES_H
#include <string>
#include <vector>
#include "all.h"
#include "Vec3.h"
#include "World.h"
#include "Field.h"
#include "Rnd.h"

class Particle{ //used os macroparticle
protected:
public:
    type_calc3 pos;  //position
    type_calc3 vel;  //velocity
    type_calc macro_weight;     //macroparticle mass

    /*constructors*/
    Particle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight) noexcept;
    Particle(type_calc3 pos, type_calc3 vel, type_calc macro_weight) noexcept;

};

class Species{
protected: 
    World&                world;      //reference to world object
    std::vector<Particle> particles;  //vector of particles

public:
    const std::string name;    //name of (micro)particle
    const type_calc   mass;    //mass of single (micro)particle
    const type_calc   charge;  //charge of sinlge (micro)particle
    const type_calc   mpw0;    //deafult macroparticle weight

    Field<type_calc> den;  //number density;
    Field<type_calc> den_avg; //average number density;

    /*constructors*/
    Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0);

    /*methods*/
    size_t getNumParticles() const;
    void advance(Species& neutrals, Species& spherium);
    void computeNumberDensity();
    void addParticle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight);
    void addParticle(type_calc3 pos, type_calc3 vel, type_calc macro_weight);
    void addParticle(type_calc3 pos, type_calc3 vel);
    //void addParticle(Particle particle); //not used, to do?
    void loadParticleBox(type_calc3 x1, type_calc3 x2, type_calc num_den, type_calc num_macro);
    void loadParticleBoxQS(type_calc3 x1, type_calc3 x2, type_calc num_den, int3 num_macro);
    type_calc  getMicroCount();   // calculate and return all number of micro particles
    type_calc3 getMomentum();     // calculate and return momentum
    type_calc  getKE();           // calculate and return kinetic energy
    type_calc  getPE();           // calculate and return potential energy
    void       updateAverages();  // update all average values

    //for advance
    type_calc3 sampleReflectedVelocity(const type_calc3& pos, const type_calc len, const type_calc3& n) const;
    type_calc sampleVth(const type_calc T) const;
};

extern Rnd rnd;
#endif