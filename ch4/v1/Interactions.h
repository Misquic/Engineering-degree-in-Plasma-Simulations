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

#endif