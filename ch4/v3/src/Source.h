#ifndef SOURCE_H
#define SOURCE_H
#include <functional>
#include "Species.h"
#include "World.h"
#include "Vec3.h"
#include "all.h"
#include "Rnd.h"

extern Rnd rnd;
/*ideas:
    make it choose a face of inlet, but overload or sth so its not computing consuming*/
class Source{
protected:
    Species& sp;     //reference to injected species
    World&   world;  //reference to world
    const type_calc v_drift;  //mean drift velocity
    const type_calc den;      // injection density
    int inlet_face_index = -1;
    const type_calc area_frac;
    //std::function<void()> sample_strategy;

    //precalc:
    type_calc3 dx;
    type_calc3 x0;
    type_calc Lx; //length x
    type_calc Ly; //length y
    type_calc Lz; //length y
    type_calc A; //surface area
    type_calc num_micro; //number of microparticles
    
    virtual void setSampleStrategy() noexcept = 0;
public:
    Source(Species& species, World& world, type_calc v_drift, type_calc den, std::string inlet_face, type_calc area_frac = 1.0) noexcept;
    virtual ~Source() noexcept = default;
    virtual void sample() const noexcept = 0;

    
};


class ColdBeamSource: public Source{
protected:
    std::function<void()> sample_strategy;

    virtual void setSampleStrategy() noexcept override;
public:
    /*constructor*/
    ColdBeamSource(Species& species, World& world, type_calc v_drift, type_calc den, std::string inlet_Face = "-z", type_calc area_frac = 1.0) noexcept;

    /*methods*/
    virtual void sample()const noexcept override;
};

class WarmBeamSource: public Source{
protected:
    std::function<void()> sample_strategy;
    type_calc T = 0;

    virtual void setSampleStrategy() noexcept override;
public:
    /*constructor*/
    WarmBeamSource(Species& species, World& world, type_calc v_drift, type_calc den, type_calc T, std::string inlet_Face = "-z", type_calc area_frac = 1.0) noexcept;

    /*methods*/
    virtual void sample()const noexcept override;
};

#endif