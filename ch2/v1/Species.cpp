#include "Species.h"
#include "Rnd.h"

/*Particle constructors*/
Particle::Particle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight) noexcept: pos{x, y, z}, vel{u, v, w}, macro_weight{macro_weight} {
};
Particle::Particle(type_calc3 pos, type_calc3 vel, type_calc macro_weight) noexcept: pos{pos}, vel{vel}, macro_weight{macro_weight} {
};

/*Species constructors*/
Species::Species(std::string name, type_calc mass, type_calc charge, World& world): name{name}, mass{mass}, charge{charge}, world{world}, den{world.ni, world.nj, world.nk} {
};

/*Species methods*/
size_t Species::getNumParticles() const{
    return particles.size();
};
void Species::advance(){
    type_calc dt = world.getDt();
    type_calc3 x0 = world.getX0();
    type_calc3 xm = world.getXm();

    type_calc3 lc{};
    type_calc3 ef_part{};
    for(Particle& part : particles){
        lc = world.XtoL(part.pos);
        
        ef_part = world.ef.gather(lc);

        part.vel += ef_part*(dt*charge/mass);

        // if(ef_part[0]*(dt*charge/mass) > 0.0001){
        //     std::cout << part.pos << " ";
        //     std::cout << part.vel*dt << " ";
        // }

        part.pos += part.vel*dt;

        // if(ef_part[0]*(dt*charge/mass) > 0.0001){
        //     std::cout << part.pos << " \n";
        // }

        //boundary check, reflective
        for(int i = 0; i < 3; i++){
            if(part.pos[i]<x0[i]){
                part.pos[i] = 2.0 * x0[i] - part.pos[i];
                part.vel[i]*= -1.0;
            }
            else if(part.pos[i]>=xm[i]){
                part.pos[i] = 2.0 * xm[i] - part.pos[i];
                part.vel[i]*= -1.0;
            }
        }
    }
};
void Species::computeNumberDensity(){

    den.clear();
    for(Particle& part : particles)
    {
        type_calc3 lc = world.XtoL(part.pos);
        den.scatter(lc, part.macro_weight);
    }

    den/=world.node_vol;
};
void Species::addParticle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight){
    addParticle({x, y, z}, {u, v, w}, macro_weight);
};
void Species::addParticle(type_calc3 pos, type_calc3 vel, type_calc macro_weight){
    if (!world.inBounds(pos)) return;

    type_calc3 lc = world.XtoL(pos);

    type_calc3 ef_part = world.ef.gather(lc);

    vel -= charge/mass*ef_part*(0.5 * world.getDt());

    particles.emplace_back(pos, vel, macro_weight);
};
void Species::loadParticleBox(type_calc3 x1, type_calc3 x2, type_calc num_den_, type_calc num_macro){
    /*x1 - starting vertex of box
      x2 - opposite vertex of box
      num den - 
      num_macro_particles - number of macro particles*/
    type_calc box_vol = (x2[0]-x1[0])*(x2[1]-x1[1])*(x2[2]-x1[2]);
    type_calc num_micro = num_den_ * box_vol;
    type_calc macro_weight = num_micro/num_macro; 

    std::cout << "Loading particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weigh of macroparticle: " << macro_weight << "\n";

    particles.reserve(num_macro);

    extern Rnd rnd;
    type_calc3 pos;
    type_calc3 vel;
    for(int i = 0; i < num_macro; i++){
        pos = {rnd(x1[0], x2[0]), rnd(x1[1], x2[1]), rnd(x1[2], x2[2])};
        vel = {0, 0, 0}; //stationary particle;
        // if(i > 50 && i < 55){
        //     std::cout << pos << "\n";
        // }

        addParticle(pos, vel, macro_weight);
    }
};
void Species::loadParticleBoxQS(type_calc3 x1, type_calc3 x2, type_calc num_den_, int3 num_macro){ //QS - quiet start
    /*x1 - starting vertex of box
      x2 - opposite vertex of box
      num den - 
      num_macro - number of macro particles per axis*/
    type_calc box_vol = (x2[0]-x1[0])*(x2[1]-x1[1])*(x2[2]-x1[2]);
    int num_macro_tot = (num_macro[0]-1)*(num_macro[1]-1)*(num_macro[2]-1);
    type_calc num_micro = num_den_ * box_vol;
    type_calc macro_weight = num_micro/num_macro_tot; 

    std::cout << "Quiet loading (uniform spacing) particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << "\nwith number of macroparticles: " << num_macro_tot << " with grid of macroparticles: " << num_macro << " and mass of macroparticle: " << macro_weight << "\n\n";

    particles.reserve(num_macro_tot);

    //particle grid spacing
    type_calc di = (x2[0] - x1[0])/(num_macro[0]-1);
    type_calc dj = (x2[1] - x1[1])/(num_macro[1]-1);
    type_calc dk = (x2[2] - x1[2])/(num_macro[2]-1);

    type_calc3 pos{};
    type_calc3 vel{};
    type_calc w{}; //number weight for mass of macroparticle to counterpart not equal cell volumes

    for(int i = 0; i < num_macro[0]; i++){
        for(int j = 0; j < num_macro[1]; j++){
            for(int k = 0; k < num_macro[2]; k++){

                pos = {x1[0] + i*di,
                       x1[1] + j*dj, 
                       x1[2] + k*dk};
                vel = {0, 0, 0}; //stationary particle;

                //shifting of particles that are on max faces of grid back to domain
                if(pos[0] == x2[0]) pos[0] -= 1e-4*di;
                if(pos[1] == x2[1]) pos[1] -= 1e-4*dj;
                if(pos[2] == x2[2]) pos[2] -= 1e-4*dk;

                w = 1;
                if(i == 0 || i == num_macro[0]-1) w*=0.5;
                if(j == 0 || j == num_macro[1]-1) w*=0.5;
                if(k == 0 || k == num_macro[2]-1) w*=0.5; 

                addParticle(pos, vel, macro_weight*w);
            }
        }
    }
};
type_calc Species::getMicroCount(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc micro_count{};
    for(Particle& part : particles){
        micro_count += part.macro_weight;
    }
    return micro_count;
};
type_calc3 Species::getMomentum(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc3 mom{};
    for(Particle& part : particles){
        mom += part.macro_weight * part.vel;
    }
    return mom * mass;
};
type_calc Species::getKE(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc ke{};
    for(Particle& part : particles){
        ke += part.macro_weight * (part.vel * part.vel); //using vector multiplication
    }
    return ke * mass * 0.5;
};
type_calc Species::getPE(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc pe{};
    type_calc phi_part{};
    type_calc3 lc{};
    for(Particle& part : particles){
        lc = world.XtoL(part.pos);
        phi_part = world.phi.gather(lc);
        pe += part.macro_weight * phi_part;
    }
    return pe * charge;
};
