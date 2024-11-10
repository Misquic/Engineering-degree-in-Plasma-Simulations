#include "Source.h"
#include "Rnd.h"

ColdBeamSource::ColdBeamSource(Species& species, World& world, type_calc v_drift, type_calc den, std::string inlet_Face) noexcept:
 sp{species}, world{world}, v_drift{v_drift}, den{den}, inlet_Face{inlet_Face}{
    dx = world.getDx();
    x0 = world.getX0();
    Lx = dx[0] * (world.ni-1);
    Ly = dx[1] * (world.nj-1);
    A = Lx*Ly;
    num_micro = den*v_drift*A*world.getDt(); //N = n*v*A*dt
};


void ColdBeamSource::sample(){
    extern Rnd rnd;
    //num_micro = den*v_drift*A*world.getDt(); //N = n*v*A*dt 
    int num_macro = (int)(num_micro/sp.mpw0 + rnd());

    type_calc3 pos{}, vel{};
    
    for(int i = 0; i < num_macro; i++){ //Z- face
        pos = {x0[0] + rnd()*Lx, x0[1] + rnd()*Ly, x0[2]};
        vel = {0, 0, v_drift};
        sp.addParticle(pos, vel, sp.mpw0);
    }

}