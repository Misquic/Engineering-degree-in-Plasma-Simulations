#include <iostream>
#include <memory>
#include "funkc.h"
#include "Instantiator.h"
#include "all.h"
#include "World.h"
#include "Field.h"
#include "Vec3.h"
#include "Outputs.h"
#include "PotentialSolver.h"
#include "Species.h"
#include "Object.h"
#include "Source.h"
#include "Rnd.h"
#include "Interactions.h"

/*TODO: add option to turn off spherium, improve ionisation, add space varying epsilon, add instantiator back, add loading with velocity sampling*/
int main(int argc, char* argv[] ){

    #ifdef DEBUG
    std::cout << "DEBUG mode\n";
    #else
    std::cout << "RELEASE mode\n";
    #endif

    std::vector<std::string> args(argv + 1, argv+argc); // passing pointers to argv values
    if(parseArgument(args, "--help") || parseArgument(args, "--h")){
        print_help();
        return 0;
    }
    
    bool SUBCYCLING = parseArgument(args, "--subcycling", true);
    Output::modes OUTPUT = parseArgument(args, "--output", Output::modes::fields); // 0 none // 1 all // 2 only screen on // 3 only screen and fields on // 4 only screen, fields and particles on etc.

    std::vector<Species> species; //species container
    std::vector<std::unique_ptr<Source>> sources; //sources container
    std::unique_ptr<PotentialSolver> solver_ptr;
    std::unique_ptr<World> world_ptr;
    std::vector<std::unique_ptr<Interaction>> interactions;

    if(parseArgument(args, "--i")){
        /*doesn't support Objects yet
          TODO trash it? command line arguments are usable and better, but Instatniator might be better for future UI?*/
        Instantiator instantiator(world_ptr, species, solver_ptr);

    }
    else{
        std::cout << "Running without Instatniator, you can use program arguments: --subcycling [<bool>] --s_type [PCG, GS] --s_max_it [<int>] --phi [<double>] --num_ts [<int>] --dt [<double>] \n\n";

        /*Read arguments*/
        SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
        int solver_max_it = parseArgument(args, "--s_max_it", 800); //or 2e4;
        type_calc phi = parseArgument(args, "--phi", -4000.0);
        int num_ts = parseArgument(args, "--num_ts", (int)4000);
        type_calc dt = parseArgument(args, "--dt", 5e-13);
        
        /*Instantiate World*/ 
        world_ptr = std::make_unique<World>(31, 31, 61, type_calc3{-0.002, -0.002, 0.0}, type_calc3{0.002, 0.002, 0.005}); //0.4cm x 0.4cm x 0.5 cm
        world_ptr->setTimeStart(); //optional, use it where you want your start of Wall time
        world_ptr->setTime(dt, num_ts);

        /*Instantiate objects*/ 
        world_ptr->addObject<Rectangle>(type_calc3(world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getX0()[2]),
         phi, type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], world_ptr->getL()[2] * 0.1));
        
        world_ptr->addObject<Rectangle>(type_calc3(world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getXm()[2]),
         -phi, type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], world_ptr->getL()[2] * 0.1));
        
        world_ptr->computeObjectID();


        type_calc E_ion_O = 1313.9*1000/Const::N_a; //in Joules
        /*Instantiate species*/ 
        species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 5e10, E_ion_O);
        species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 2e2);  // ions //last is mpw0
        // species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e0);
        species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr, 2e2);      // electrons

        const type_calc num_den_neutrals = 1e25; // final should be 1e25 
        const type_calc num_den_ions = 1e10; //mean ion density
        const type_calc num_den_electrons = 1e19; //mean electron density
        type_calc T = 300;
        type_calc T_eles = 10*T;

        // type_calc side_z = (world_ptr->getXm()[2] - world_ptr->getX0()[2])*0.9;
        // ZROBIC dodać więcej czynników do samplowania prędkości i sprawdzić czy to pomaga
        species[2].loadParticleBoxThermal(type_calc3{world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getXc()[2] - 0.5*world_ptr->getL()[2]*0.85}, type_calc3(world_ptr->getL()[0]*0.02, world_ptr->getL()[1]*0.02, world_ptr->getL()[2]*0.02), num_den_electrons, T_eles);

        species[0].loadParticleBoxThermal(world_ptr->getXc(), type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], world_ptr->getL()[2]*0.9), num_den_neutrals, T);
        
        /* Instantiate interactions*/
        interactions.emplace_back(std::make_unique<MC_MEX_IONIZATION>(species[0], species[1], species[2], *world_ptr));
        interactions.emplace_back(std::make_unique<DSMC_MEX>(species[2], *world_ptr));

        //interactions.emplace_back(std::make_unique<DSMC_MEX_IONIZATION>(species[1], species[2], species[3], *world_ptr));

        // type_calc rate_0_1 = 1e-2;
        // type_calc rate_1_2 = 1e-3;
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[0], species[1], species[3], *world_ptr, rate_0_1));
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[1], species[2], species[3], *world_ptr, rate_1_2));

        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));


        /*Instantiate solver*/ 
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        solver_ptr->setReferenceValues(0, num_den_ions, T);

        std::cout << "\nOptions:\n";
        std::cout << "Sybcycling: " << SUBCYCLING << "\n";
        std::cout << "output mode: " << OUTPUT << "\n";
        std::cout << "Solver Type: " << solver_type << "\n";
	    std::cout << "GS Solver max iterations: " << solver_ptr->get_GS_max_it() <<"\n";
        if(solver_type==PCG){
    	    std::cout << "PCG Solver max iterations: " << solver_ptr->get_PCG_max_it() <<"\n";
        }
	    std::cout << "Potential at one electrode: " << phi << " V\n" << std::endl;

    }
    Species& neutral_oxygen = species[0];
    Species& one_oxygen = species[1];
    // Species& two_oxygen = species[2];
    Species& electrons = species[2];


    //////////////
    
    for (Species &sp : species) {
        std::cout << sp.name << " has " << sp.getNumParticles() << " particles" << std::endl;
        sp.computeMacroParticlesCount();
    }

    //Output::fields(*world_ptr, species);
    solver_ptr->solve();
    solver_ptr->computeEF();
    Output::fieldsOutput(*world_ptr, species);

    //world_ptr->setTimeStart();
    type_calc dt = world_ptr->getDt();
    type_calc start_time_meas = 0;
    type_calc time_interactions = 0;
    type_calc time_advances = 0;
    type_calc time_computes = 0;
    while(world_ptr->advanceTime()){

        for(std::unique_ptr<Source>& source: sources){
            source->sample();
        }

        for(std::unique_ptr<Interaction>& interaction: interactions){
            interaction->apply(world_ptr->getDt());
        }

        if(SUBCYCLING == true){
            for(Species&  sp: species){//TODO correct it so for species neutrals refer to species neutrals
                if(sp.name == electrons.name){ //implemented subcycling //electrons
                    sp.advanceElectrons(dt);
                    sp.computeNumberDensity();
                    sp.sampleMoments();
                    sp.computeMacroParticlesCount();
                
                    // std::cout << "moved electrons ";
                }
                else if(sp.charge!= 0 && world_ptr->getTs()%10 == 0){ //ions
                    // sp.advance(neutral_oxygen, neutral_oxygen, dt*10); //with sputtering material
                    sp.advanceNoSputtering(neutral_oxygen, dt*10); //without sputtering material
                    sp.computeNumberDensity();
                    sp.sampleMoments();
                    sp.computeMacroParticlesCount();
                    // std::cout << " moved ions ";
                }
                else if(world_ptr->getTs()%100 == 0){ // neutrals 
                    //sp.advanceNoSputtering(neutral_oxygen, dt*50); //without sputtering material, neutrals can't do it
                    sp.computeNumberDensity();
                    sp.sampleMoments();
                    sp.computeMacroParticlesCount();
                    // std::cout << " moved neutrals ";
                }
                // sp.advance(neutral_oxygen, neutral_oxygen, dt); //normal
            }
        }
        else{
            for(Species&  sp: species){//TODO correct it so for species neutrals refer to species neutrals
                if(sp.name == electrons.name){ //implemented subcycling //electrons
                    sp.advanceElectrons(dt);
                }
                else{ //ions
                    // sp.advanceSputtering(neutral_oxygen, neutral_oxygen, dt); //with sputtering material
                    sp.advanceNoSputtering(neutral_oxygen, dt); //no sputtering
                }
                sp.computeNumberDensity();
                sp.sampleMoments();
                sp.computeMacroParticlesCount();
            }
        }

        if(world_ptr->getTs()>5){ //moved before next section
            for(Species&  sp: species){
                sp.updateAverages();
            }
        }

        world_ptr->computeChargeDensity(species);
        solver_ptr->solve();
        solver_ptr->computeEF();


        if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::diagnostics ){
            Output::diagOutput(*world_ptr, species);
        }
        if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::screen ){
            Output::screenOutput(*world_ptr, species);
        }
        int ts = world_ptr->getTs();
        if(ts%10 == 0 || world_ptr->isLastTimeStep()){ //|| (ts > 140 && ts < 160)){
            if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::fields ){
                Output::fieldsOutput(*world_ptr, species);
            }
            if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::particles ){
                Output::particlesOutput(*world_ptr, species, 1000);
            }
            if(world_ptr->getWallTime()>60){
                std::cout << "Time taken so far " << world_ptr->getWallTime()/60 << " minutes." << std::endl;
            }else{
                std::cout << "Time taken so far " << world_ptr->getWallTime() << " seconds." << std::endl;
    }
        }
        std::cout << std::flush <<"\n";

       
    }
    if(world_ptr->getWallTime()>60){
        std::cout << "Simulation took " << world_ptr->getWallTime()/60 << " minutes." << std::endl;
    }else{
        std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
    }
    return 0;
}

