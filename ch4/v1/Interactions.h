#ifndef INTERACTIONS_H
#define INTERACTIONS_H
#include "all.h"
#include "Species.h"
#include "World.h"

class Interaction{
public:
    virtual void apply(type_calc dt) noexcept = 0;
    virtual ~Interaction() noexcept = default;
};

class ChemistryIonize: public Interaction{
protected:
    Species& neutrals;
    Species& ions;
    World& world;
    type_calc rate;
    
public:
    ChemistryIonize(Species& neutrals, Species& ions, World& world, type_calc rate) noexcept;
    virtual void apply(type_calc dt) noexcept override;
    Field<type_calc> dNi_field;
};

class DSMC_MEX: public Interaction{
protected:
    Species& species;
    World& world;

    type_calc m_reduced;
    type_calc c[4]; //for sigma evaluation according to Bird's

	type_calc sigma_v_rel_max = 1e-14;	//some initial value
    type_calc evaluateSigma(type_calc g_rel);
    void collide(type_calc3& vel1, type_calc3& vel2, type_calc mass1, type_calc mass2);
public:
    DSMC_MEX(Species& species, World& world) noexcept;
    virtual void apply(type_calc dt) noexcept override;

};

#endif