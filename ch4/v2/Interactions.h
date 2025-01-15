#ifndef INTERACTIONS_H
#define INTERACTIONS_H
#include "all.h"
#include "Species.h"
#include "World.h"
#include <functional>
#include <map>

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

class ChemistryIonizeElectrons: public Interaction{ //clas to use when directly simulating electrons
protected:
    Species& neutrals;
    Species& ions;
    Species& electrons;
    World& world;
    type_calc rate;
    
public:
    ChemistryIonizeElectrons(Species& neutrals, Species& ions, Species& electrons, World& world, type_calc rate) noexcept;
    virtual void apply(type_calc dt) noexcept override;
    Field<type_calc> dNi_field;
};

class DSMC_MEX: public Interaction{
protected:
    Species& species1;
    Species& species2;
    World& world;
    std::function<void(type_calc)> apply_strategy;

    //precalculate
    type_calc m_reduced;
    type_calc mass1;
    type_calc mass2;
    type_calc sum_mass;
    type_calc dv;

    type_calc c[4]; //for sigma evaluation according to Bird's

	type_calc sigma_v_rel_max = 1e-14;	//some initial value
    type_calc evaluateSigma(type_calc v_rel);
    void collide(type_calc3& vel1, type_calc3& vel2);
public:
    DSMC_MEX(Species& species, World& world) noexcept;
    DSMC_MEX(Species& species1, Species& species2, World& world);
    virtual void apply(type_calc dt) noexcept override;
    virtual void applyOneSpecies(type_calc dt) noexcept;
    virtual void applyTwoSpecies(type_calc dt) noexcept;
};

class DSMC_MEX_Ionization: public Interaction{
protected:
    Species& neutrals; //just called so, can be charged, but one less than ions
    Species& ions;
    Species& electrons;
    World& world;
    std::function<void(type_calc)> apply_strategy;

    type_calc E_ion_times_10; //in Joules

    //precalculate
    type_calc m_reduced;
    type_calc mass_neus;
    type_calc mass_eles;
    type_calc sum_mass;
    type_calc dv;
    int num_cells;

    type_calc c[4]; //for sigma evaluation according to Bird's

	// type_calc sigma_v_rel_max = 1e-14;	//some initial value
	type_calc sigma_v_rel_max = 1e-14;	//some initial value
    type_calc evaluateSigma(type_calc v_rel);
    void collide(type_calc3& vel1, type_calc3& vel2, bool& ionize, type_calc& E_new);
public:
    DSMC_MEX_Ionization(Species& neutrals, Species& ions, Species& electrons, World& world);

    virtual void apply(type_calc dt) noexcept override;
};

class MC_MEX_Ionization: public Interaction{
protected:
    Species& neutrals; //just called so, can be charged, but one less than ions
    Species& ions;
    Species& electrons;
    World& world;
    std::function<void(type_calc)> apply_strategy;

    type_calc E_ion_times_10; //in Joules
    type_calc E_ion_J; //in Joules
    type_calc E_ion_eV; //in eV

    //precalculate
    type_calc m_reduced;
    type_calc mass_neus;
    type_calc mass_eles;
    type_calc sum_mass;
    type_calc dv;
    type_calc inv_dv;
    type_calc dt;
    int num_cells;
    type_calc E_precalc_rel_J;
    type_calc E_precalc_rel_eV;
    type_calc E_precalc_ele_eV = Const::m_e *0.5/Const::q_e;
    type_calc two_qe_to_me = 2*Const::q_e/Const::m_e;

    type_calc c[4]; //for sigma evaluation according to Bird's
    type_calc B_inc;

	// type_calc sigma_v_rel_max = 1e-14;	//some initial value
	type_calc sigma_v_rel_max = 1e-14;	//some initial value
    type_calc evaluateSigmaColl(type_calc v_rel);
    type_calc evaluateSigmaIon(type_calc v_rel);
    std::map<type_calc, type_calc> data_sigma;
    void collide(type_calc3& vel_atom, type_calc3& vel_ele, bool& ionize, type_calc3& vel_new, int N_atoms_macro, type_calc sigma_coll);
    void newVelocityElecton(type_calc3& vel, type_calc E, type_calc3& vel_inc_unit); //IONIZE_1
    void newVelocityElecton(type_calc3& vel, type_calc E); //IONIZE_2
public:
    MC_MEX_Ionization(Species& neutrals, Species& ions, Species& electrons, World& world, std::string coll_data_path = "data/Oxygen_momentum_transfer.txt");

    virtual void apply(type_calc dt) noexcept override;
};

#endif