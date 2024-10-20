#ifndef WORLD_H
#define WORLD_H
#include "all.h"
#include "Field.h"
#include "Vec3.h"
#include <chrono>


class Species;

class World{
protected:
    type_calc3 x0;  //mesh origin
    type_calc3 dx;  //cell spacing
    type_calc3 xm;  //mesh max bound
    type_calc3 xc;  //mesh center
    type_calc dt = 1e-4;     //time step
    int num_ts = 0;     //time step
    int ts = -1;     //time step
    type_calc time = 0;
    const std::chrono::time_point<std::chrono::high_resolution_clock> time_start;

    /*protected methods*/
    void computeNodeVolumes();

public:
    const int nn[3];       //number of nodes
    const int ni, nj, nk;  //number of nodes in individual directions

    //TO DO take this to protected and add friends which use that Fields
    Field<type_calc>  phi;       //potential
    Field<type_calc>  rho;       //charge density
    Field<type_calc>  node_vol;  //node volumes
    Field<type_calc3> ef;        //electric field
    
    /*constructors*/
    World(int ni, int nj, int nk);
    World(int ni, int nj, int nk, type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2);
    World(int ni, int nj, int nk, type_calc3 vec1, type_calc3 vec2);

    /*methods*/
    void setExtents(type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2);
    void setExtents(type_calc3 vec1, type_calc3 vec2);
    type_calc3 getX0() const;
    type_calc3 getDx() const;
    type_calc3 getXm() const;
    type_calc3 getXc() const;
    type_calc3 XtoL(type_calc3 x) const;
    type_calc getPE() const;

    //type_calc3 getTs() const;
    void computeChargeDensity(std::vector<Species>& species);
    bool inBounds(type_calc3 pos);

    //time
    void setTime(type_calc dt, int num_ts);
    type_calc getDt() const;
    int getTs() const;
    type_calc getTime() const;
    bool isLastTimeStep() const;
    bool advanceTime();
    type_calc getWallTime();


};


#endif