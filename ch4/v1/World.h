#ifndef WORLD_H
#define WORLD_H
#include <chrono>
#include <memory>
#include <tuple>
#include "all.h"
#include "Field.h"
#include "Vec3.h"
#include "Object.h"



class Species;

class World{
protected:
    type_calc3 x0;  //mesh origin
    type_calc3 dx;  //cell spacing
    type_calc3 xm;  //mesh max bound
    type_calc3 xc;  //mesh center

    /*time*/
    type_calc dt = 1e-4;     //time step
    int num_ts = 0;     //time step
    int ts = -1;     //time step
    type_calc time = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> time_start;

    /*steady state*/
    bool steady_state = false;
    type_calc last_mass = 0;
    type_calc last_mom = 0;
    type_calc last_en = 0;

    
    std::vector<std::shared_ptr<Object>> objects;

    /*protected methods*/
    void computeNodeVolumes();

public:
    const int nn[3];       //number of nodes
    const int ni, nj, nk;  //number of nodes in individual directions
    const int nv;          //volume in nodes

    //TODO take this to protected and add friends which use that Fields
    Field<type_calc>  phi;        //potential
    Field<type_calc>  rho;        //charge density
    Field<type_calc>  node_vol;   //node volumes
    Field<type_calc3> ef;         //electric field
    Field<int>        object_id;  //object id flag to flag fixed nodes //calculating that via Objects not World (unlike the book)
    Field<type_calc>  object_phi;  //object phi for fixed nodes //calculating that via Objects not World (unlike the book)
    Field<int>  node_type;  //node_type
    
    /*constructors*/
    World(int ni, int nj, int nk);
    World(int3 nn);
    World(int ni, int nj, int nk, type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2);
    World(int ni, int nj, int nk, type_calc3 vec1, type_calc3 vec2);
    World(const World& other) noexcept;
    World(World&& other) noexcept;

    /*methods*/
    void setExtents(type_calc x1, type_calc y1, type_calc z1, type_calc x2, type_calc y2, type_calc z2);
    void setExtents(type_calc3 vec1, type_calc3 vec2);
    type_calc3 getX0() const;
    type_calc3 getDx() const;
    type_calc3 getXm() const;
    type_calc3 getXc() const;
    type_calc getPE() const;
    type_calc3 XtoL(type_calc3 x) const;   //converts position to logical coordinates
    type_calc3 LtoX(type_calc3 lc) const;  //converts logical coordinates to position
    type_calc3 LtoX(int3 lc) const;  //converts logical coordinates to position
    type_calc3 LtoX(int i, int j, int k) const;  //converts logical coordinates to position
    bool steadyState(std::vector<Species>& species);

    //type_calc3 getTs() const;
    void computeChargeDensity(std::vector<Species>& species);
    bool inBounds(type_calc3 pos);
    void addInlet(std::string face = "-z");

    /*Objects*/
    void computeObjectID();
    int inObject(const type_calc3& pos) const; // returns index of object pos is in, starting from 1, 0 means pos isn't in any object!!!
    void lineIntersect(const type_calc3& x1, const type_calc3& x2, const int in_object, type_calc& tp, type_calc3& intersection_point, type_calc3& n) const; //change it so it takes references of values in stead of returning tuple?
    template<class T, class... Args>
    void addObject(Args&&... args);
    std::string printObjects()const;
    void printObjects(std::ostream& out)const;

    /*time*/
    void setTime(type_calc dt, int num_ts);
    type_calc getDt() const;
    int getTs() const;
    type_calc getTime() const;
    bool isLastTimeStep() const;
    bool advanceTime();
    type_calc getWallTime();
    void setTimeStart();

    World& operator=(const World& other);
    World& operator=(World&& other);
};

template <class T, class... Args>
void World::addObject(Args&&... args){
    static_assert(std::is_base_of<Object, T>::value, "T must derive from Object");
    try{
        objects.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error adding object, supposedly you used wrong arguments that do not match any of constructor of Object or derived class: " << e.what() << '\n';
    } catch (...) {
        std::cerr << "Unknown error occurred while adding object\n";
    }

}

#endif