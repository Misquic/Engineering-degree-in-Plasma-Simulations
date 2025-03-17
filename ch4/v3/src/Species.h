#ifndef SPECIES_H
#define SPECIES_H
#include <string>
#include <vector>
#include <unordered_map>
#include "all.h"
#include "Vec3.h"
#include "World.h"
#include "Field.h"
#include "Rnd.h"

class Particle{ //used as macroparticle
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
    std::vector<std::vector<int>> sorted_particles_indexes;  //vector of sorted indexes per cell
    bool sorted = false;

    Field<type_calc>  n_sum;                      //number of particles
    Field<type_calc3> nv_sum;                     //number * velocity
    Field<type_calc>  nuu_sum, nvv_sum, nww_sum;  // number*vel squared

public:
    const std::string name;    //name of (micro)particle
    const type_calc   mass;    //mass of single (micro)particle
    const type_calc   charge;  //charge of sinlge (micro)particle
    const type_calc   mpw0;    //deafult macroparticle weight
    const type_calc   E_ion = -666;   //energy of ionisation of single (micro)particle in Joules, nan means species has no way of ionizating

    Field<type_calc>  den;               //number density;
    Field<type_calc>  den_avg;           //average number density;
    Field<type_calc>  T;                 //temperature
    Field<type_calc3> vel;               //stream velocity
    Field<type_calc>  macro_part_count;  //field of number ofmacroparticles per cell

    /*constructors*/
    Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0);
    Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0, type_calc E_ion);

    /*methods*/
    size_t getNumParticles() const;
    type_calc advance_time_multi = -1; //for automatic switching for instance of object, first it will measure multithreading, then serial, then switching is valid
    type_calc advance_time_serial = 0;
    void advanceNonElectron(Species& neutrals, Species& spherium, type_calc dt);
private:
    void advanceSputteringSerial(Species& neutrals, Species& spherium, type_calc dt); //with sputtering material
    void advanceNoSputteringSerial(Species& neutrals, type_calc dt); //without sputtering material
public:
    void advanceElectrons(type_calc dt); //for electrons
private:
    void advanceElectronsMultithreading(unsigned int thread_id, type_calc dt, size_t index_start, size_t index_stop, std::vector<size_t>& parts_to_delete);
    void advanceElectronsSerial(type_calc dt);
public:


    void computeNumberDensity();
    void addParticle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight);
    void addParticle(type_calc3 pos, type_calc3 vel, type_calc macro_weight);
    void addParticle(type_calc3 pos, type_calc3 vel);
    //void addParticle(Particle particle); //not used, to do?
    void loadParticleBox(type_calc3 x1, type_calc3 x2, type_calc num_den, type_calc num_macro);
    void loadParticleBoxQS(type_calc3 x1, type_calc3 x2, type_calc num_den, int3 num_macro);
    void loadParticleBoxQS(type_calc3 x1, type_calc3 x2, type_calc num_den);
    void loadParticleBoxThermal(type_calc3 x0, type_calc3 sides, type_calc num_den, type_calc T);
    void loadParticleSphereThermal(type_calc3 x0, type_calc r, type_calc num_den, type_calc T);
    void loadLangumir(type_calc lambda, type_calc A, type_calc num_den_, type_calc T);
    void move(type_calc3 dx, Rectangle& rec); // move particles that are in rec dx 
    void move(type_calc3 dx);
    void move(type_calc lambda, type_calc A);
    bool setSorted(bool is_sorted) noexcept;
    bool isSorted() noexcept;

    /*for diagnostics and outputs*/
    type_calc  getMicroCount();                         // calculate and return all number of micro particles
    type_calc3 getMomentum();                           // calculate and return momentum
    type_calc  getKE();                                 // calculate and return kinetic energy
    type_calc  getPE();                                 // calculate and return potential energy
    void       updateAverages();                        // update all average values
    void       sampleMoments();                         // sample moments od velocity distribution function (n, nv, nuu, nvv, nww)
    void       computeGasProperties();                  // use sampled momentss to compute velocity and temperature
    void       clearSamples();                          // after writing to file clear Samples
    void       computeMacroParticlesCount();            // compute number of macroparticles per cell

    const std::vector<Particle>& getConstPartRef() const;  //get constant reference to particles (following encapsulation principles)
    std::vector<Particle>& getPartRef();  //get constant reference to particles (following encapsulation principles)
    const Particle& getConstPartRef(int i) const;          //get constant reference to particle (following encapsulation principles)

    //for advance
    type_calc3 sampleReflectedVelocity(const type_calc3& pos, const type_calc v_mag1, const type_calc3& n) const;
    type_calc sampleVth(const type_calc T) const;
    type_calc3 sampleV3th(const type_calc T) const; 
    void sampleVthVariableMpw(const type_calc T, type_calc& set_vel, type_calc& set_mpw) const; //TODO implemet it and particle merge?


    //for sorting
    std::vector<std::vector<Particle*>> sortPointers();
    //std::vector<std::vector<std::unique_ptr<Particle>>> sort_pointers_unique();
    std::vector<std::vector<int>> sortIndexes();
    std::vector<std::vector<int>> getSortedIndexes(); //TODO change it to ref
    void sortIndexes(std::vector<std::vector<int>>& particles_in_cell);
    std::unordered_map<int, std::vector<int>> mapIndexes();
    void mapIndexes(std::unordered_map<int, std::vector<int>>& map_particles_in_cell);
    void mapIndexes(std::unordered_map<int, std::vector<int>>& map_new_particles_in_cell, std::unordered_map<int, std::vector<int>>& already_mapped_other_species_particles);
    
    //for merge

    void merge();
    bool shouldUseMap();
protected:
    const type_calc3 m_vel_grid_n{15.0};
    std::vector<std::vector<std::vector<std::vector<int>>>> sortVelocitiesInCell(const std::vector<int>& particles_in_cell);


};

extern Rnd rnd;
#endif