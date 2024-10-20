#ifndef SPECIES_H
#define SPECIES_H
#include <string>
#include <vector>
#include "all.h"
#include "Vec3.h"
#include "World.h"
#include "Field.h"

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

    Field<type_calc> den;  //number density;

    /*constructors*/
    Species(std::string name, type_calc mass, type_calc charge, World& world);

    /*methods*/
    size_t getNumParticles() const;
    void advance();
    void computeNumberDensity();
    void addParticle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight);
    void addParticle(type_calc3 pos, type_calc3 vel, type_calc macro_weight);
    //void addParticle(Particle particle); //not used, to do?
    void loadParticleBox(type_calc3 x1, type_calc3 x2, type_calc num_den, type_calc num_macro);
    void loadParticleBoxQS(type_calc3 x1, type_calc3 x2, type_calc num_den, int3 num_macro);
    type_calc getMicroCount();
    type_calc3 getMomentum();
    type_calc getKE();
    type_calc getPE();
};
#endif