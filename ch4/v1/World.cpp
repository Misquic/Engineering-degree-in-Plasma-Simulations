#include <map>
#include <sstream>
#include "funkc.h"
#include "World.h"
#include "Species.h"

/*constructors*/
World::World(int ni, int nj, int nk): nn{ni, nj, nk}, ni{ni}, nj{nj}, nk{nk}, nv{ni*nj*nk}, phi(nn), rho(nn), node_vol(nn), ef(nn), object_id(nn), object_phi(nn), node_type(nn),
 time_start{std::chrono::high_resolution_clock::now()}{
};
World::World(int3 nn): nn{nn[0], nn[1], nn[2]}, ni{nn[0]}, nj{nn[1]}, nk{nn[2]}, nv{ni*nj*nk}, phi(nn), rho(nn), node_vol(nn), ef(nn), object_id(nn), object_phi(nn), node_type(nn),
 time_start{std::chrono::high_resolution_clock::now()}{
};
World::World(int ni, int nj, int nk, type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2): nn{ni, nj, nk},
 ni{ni}, nj{nj}, nk{nk}, nv{ni*nj*nk}, phi(nn), rho(nn), node_vol(nn), ef(nn), object_id(nn), object_phi(nn), node_type(nn), time_start{std::chrono::high_resolution_clock::now()}{
    setExtents(x1, y1, z1, x2, y2, z2);
};
World::World(int ni, int nj, int nk, type_calc3 vec1, type_calc3 vec2): nn{ni, nj, nk}, ni{ni}, nj{nj}, nk{nk}, nv{ni*nj*nk}, phi(nn), rho(nn), node_vol(nn),
 ef(nn), object_id(nn), object_phi(nn), node_type(nn), time_start{std::chrono::high_resolution_clock::now()}{
    setExtents(vec1, vec2);
};
World::World(const World& other) noexcept: nn{other.ni, other.nj, other.nk}, ni{other.ni},
 nj{other.nj}, nk{other.nk}, nv{other.nv}, phi(other.phi), rho(other.rho), node_vol(other.node_vol),
 ef(other.ef), object_id(other.object_id), object_phi(other.object_phi), node_type(other.node_type), time_start{std::chrono::high_resolution_clock::now()}, x0{other.x0}, dx{other.dx},
 xm{other.xm}, xc{other.xc}, dt{other.dt}, num_ts{other.num_ts}, ts{other.ts}, time{other.time} {
};
World::World(World&& other) noexcept:
 nn{other.ni, other.nj, other.nk},
 ni{other.ni},
 nj{other.nj},
 nk{other.nk},
 nv{other.nv},
 phi(std::move(other.phi)), 
 rho(std::move(other.rho)), 
 node_vol(std::move(other.node_vol)),
 ef(std::move(other.ef)), 
 object_id(std::move(other.object_id)), 
 object_phi(std::move(other.object_phi)), 
 node_type(std::move(other.node_type)),
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
};
type_calc3 World::XtoL(type_calc3 x) const{ //L - length from begining
    type_calc3 lc = (x - x0)/dx;

    return lc; //in world coz x0 and dx used
};
type_calc3 World::LtoX(type_calc3 lc) const{ //converts logical coordinates to position
static type_calc3 ret{};
    for(int i = 0; i < 3; i++){
        ret[i] = x0[i] + lc[i] * dx[i];
    }
    return ret;
};
type_calc3 World::LtoX(int3 lc) const{ //converts logical coordinates to position
    return LtoX(type_calc3(lc[0], lc[1], lc[2]));
};
type_calc3 World::LtoX(int i, int j, int k) const{ //converts logical coordinates to position
    return LtoX(type_calc3(i, j, k));
};
bool World::steadyState(std::vector<Species>& species){
    if(steady_state) return true;

    type_calc tot_mass = 0;
    type_calc tot_mom = 0;
    type_calc tot_en = getPE();
    for( Species& sp: species){
        tot_mass += sp.getMicroCount();
        tot_mom += sp.getMomentum().length();
        tot_en += sp.getKE();
    }
    const type_calc tol = 1e-3;
    if(std::fabs((tot_mass - last_mass)/tot_mass) < tol && std::fabs(tot_mom - last_mom)/tot_mom < tol && std::fabs(tot_en - last_en)/tot_en < tol){
        steady_state = true;
        std::cout << "\nSteady state reached at time step " << ts << std::endl;
    }
    last_mass = tot_mass;
    last_mom = tot_mom;
    last_en = tot_en;
    return steady_state;
};


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
void World::addInlet(std::string face){
    static std::map<std::string, int> faces{ //map to map faces to ints to work easier later
        {"x-", 0},
        {"x+", 1},
        {"y-", 2},
        {"y+", 3},
        {"z-", 4},
        {"z+", 5},
        {"-x", 0},
        {"+x", 1},
        {"-y", 2},
        {"+y", 3},
        {"-z", 4},
        {"+z", 5},
    };
    lower(face);
    if(faces.find(face) == faces.end()){ //checking if passed wrong argument
        throw std::invalid_argument("Passed wrong face name.");
        return;
    }
    int find = faces[face];
    switch(find){
    case 0:
        for(int j = 0; j < nj; j++){
            for(int k = 0; k < nk; k++){
                object_id[0][j][k] = 2;
                phi[0][j][k] = 0;
            }
        }
        break;
    case 1:
        for(int j = 0; j < nj; j++){
            for(int k = 0; k < nk; k++){
                object_id[ni-1][j][k] = 2;
                phi[ni-1][j][k] = 0;
            }
        }
        break;
    case 2:
        for(int i = 0; i < ni; i++){
            for(int k = 0; k < nk; k++){
                object_id[i][0][k] = 2;
                phi[i][0][k] = 0;
            }
        }
        break;
    case 3:
        for(int i = 0; i < ni; i++){
            for(int k = 0; k < nk; k++){
                object_id[i][nj-1][k] = 2;
                phi[i][nj-1][k] = 0;
            }
        }
        break;
    case 4:
        for(int i = 0; i < ni; i++){
            for(int j = 0; j < nj; j++){
                object_id[i][j][0] = 2;
                phi[i][j][0] = 0;
            }
        }
        break;
    case 5:
        for(int i = 0; i < ni; i++){
            for(int j = 0; j < nj; j++){
                object_id[i][j][nk-1] = 2;
                phi[i][j][nk-1] = 0;
            }
        }
        break;
    default:
        throw std::invalid_argument("Something went wrong, addInlet default for switch");
        break;
    }
};
void World::computeObjectID(){
    for(const std::shared_ptr<Object>& obj_ptr: objects){
        for(int i = 0; i < this->ni; i++){
            for(int j = 0; j < this->nj; j++){
                for (int k = 0; k < this->nk; k++){
                    if(obj_ptr->inObject(this->LtoX(i, j, k))){
                        object_id[i][j][k] = 1;
                        //object_phi[i][j][k] = obj_ptr->getPhi();
                        phi[i][j][k] = obj_ptr->getPhi();
                    }
                }
            }
        }
    }
};
bool World::inObject(const type_calc3& pos) const{
    for(const std::shared_ptr<Object>& obj_ptr: objects){
        if(obj_ptr->inObject(pos)) return true;
    }
    return false;
};
std::string World::printObjects()const{
    std::stringstream ss;
    for(const std::shared_ptr<Object>& obj_ptr: objects){
        ss << *obj_ptr << "\n";
    }
    return ss.str();
};
void World::printObjects(std::ostream& out)const{
    for(const std::shared_ptr<Object>& obj_ptr: objects){
        out << *obj_ptr << "\n";
    }
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
