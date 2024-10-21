#include "World.h"
#include "Species.h"

//TODO add instantiating Field by int3 and int[3]
/*constructors*/
World::World(int ni, int nj, int nk): nn{ni, nj, nk}, ni{ni}, nj{nj}, nk{nk}, phi(nn), rho(nn), node_vol(nn), ef(nn), object_id(nn),
 time_start{std::chrono::high_resolution_clock::now()}{
};
World::World(int3 nn): nn{nn[0], nn[1], nn[2]}, ni{nn[0]}, nj{nn[1]}, nk{nn[2]}, phi(nn), rho(nn), node_vol(nn), ef(nn), object_id(nn),
 time_start{std::chrono::high_resolution_clock::now()}{
};
World::World(int ni, int nj, int nk, type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2): nn{ni, nj, nk},
 ni{ni}, nj{nj}, nk{nk}, phi(nn), rho(nn), node_vol(nn), ef(nn), object_id(nn), time_start{std::chrono::high_resolution_clock::now()}{
    setExtents(x1, y1, z1, x2, y2, z2);
};
World::World(int ni, int nj, int nk, type_calc3 vec1, type_calc3 vec2): nn{ni, nj, nk}, ni{ni}, nj{nj}, nk{nk}, phi(nn), rho(nn), node_vol(nn),
 ef(nn), object_id(nn), time_start{std::chrono::high_resolution_clock::now()}{
    setExtents(vec1, vec2);
};
World::World(const World& other) noexcept: nn{other.ni, other.nj, other.nk}, ni{other.ni},
 nj{other.nj}, nk{other.nk}, phi(other.phi), rho(other.rho), node_vol(other.node_vol),
 ef(other.ef), object_id(other.object_id), time_start{std::chrono::high_resolution_clock::now()}, x0{other.x0}, dx{other.dx},
 xm{other.xm}, xc{other.xc}, dt{other.dt}, num_ts{other.num_ts}, ts{other.ts}, time{other.time} {
};
World::World(World&& other) noexcept:
 nn{other.ni, other.nj, other.nk},
 ni{other.ni},
 nj{other.nj},
 nk{other.nk}, 
 phi(std::move(other.phi)), 
 rho(std::move(other.rho)), 
 node_vol(std::move(other.node_vol)),
 ef(std::move(other.ef)), 
 object_id(std::move(other.object_id)), 
 time_start{std::chrono::high_resolution_clock::now()}, 
 x0{other.x0}, 
 dx{other.dx},
 xm{other.xm}, 
 xc{other.xc}, 
 dt{other.dt}, 
 num_ts{other.num_ts},
 ts{other.ts}, 
 time{other.time} {
};

/*methods*/
void World::setExtents(type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2){
    x0[0] = x1;
    x0[1] = y1;
    x0[2] = z1;

    xm[0] = x2;
    xm[1] = y2;
    xm[2] = z2;

    for(int i = 0; i < 3; i++){
        dx[i] = (xm[i] -x0[i])/(nn[i] - 1); 
        xc[i] = 0.5 * (xm[i] + x0[i]);
    }
    computeNodeVolumes();
};
void World::setExtents(type_calc3 vec1, type_calc3 vec2){
    setExtents(vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
};
type_calc3 World::getX0() const{
    return type_calc3(x0);
};
type_calc3 World::getDx() const{
    return type_calc3(dx);
};
type_calc3 World::getXm() const{
    return type_calc3(xm);
};
type_calc3 World::getXc() const{
    return type_calc3(xc);
};
type_calc3 World::XtoL(type_calc3 x) const{ //L - length from begining
    type_calc3 lc{};
    lc[0] = (x[0] - x0[0])/dx[0];
    lc[1] = (x[1] - x0[1])/dx[1];
    lc[2] = (x[2] - x0[2])/dx[2];

    return lc; //in world coz x0 and dx used
};
type_calc World::getPE() const{
    type_calc pe{};
    for(int i = 0; i < ni; i++){
        for(int j = 0; j < nj; j++){
            for(int k = 0; k < nk; k++){
                pe += ef[i][j][k] * ef[i][j][k] * node_vol[i][j][k]; 
            }
        }
    }
    return 0.5*Const::eps_0*pe;
}


void World::computeChargeDensity(std::vector<Species>& species){
    rho = 0;
    for(Species& sp : species){
        //sp.computeNumberDensity();
        if(sp.charge == 0) continue;
        rho += sp.charge * sp.den;
    }
};
bool World::inBounds(type_calc3 pos){
    for (int i=0;i<3;i++)
        if (pos[i]<x0[i] || pos[i]>=xm[i]) return false;
    return true;
};

//time
void World::setTime(type_calc dt, int num_ts){
    this->dt = dt;
    this->num_ts = num_ts;
};
type_calc World::getDt() const{
    return dt;
};
int World::getTs() const{
    return ts;
};
type_calc World::getTime() const{
    return time;
};
bool World::isLastTimeStep() const{
    return ts == num_ts-1;
};
bool World::advanceTime(){
    time+=dt;
    ts++;
    return ts<=num_ts;
};
type_calc World::getWallTime(){
    auto time_end = std::chrono::high_resolution_clock::now();
    type_calc time_taken = (time_end - time_start).count()*1e-9;
    return time_taken;
};
void World::setTimeStart(){
    time_start = std::chrono::high_resolution_clock::now();
}


/*protected methods*/
void World::computeNodeVolumes(){
    type_calc base_volume = dx[0]*dx[1]*dx[2];
    type_calc calc_volume = base_volume;
    for(int i = 0; i < ni; i++){ //?write explicite without if? for performance
        for(int j = 0; j < nj; j++){
            for(int k = 0; k < nk; k++){
                calc_volume = base_volume;
                if(i == 0 || i == ni-1) calc_volume*=0.5;
                if(j == 0 || j == nj-1) calc_volume*=0.5;
                if(k == 0 || k == nk-1) calc_volume*=0.5;
                node_vol[i][j][k] = calc_volume;
            }
        }
    }
}
