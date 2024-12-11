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

    std::vector<std::string> args(argv + 1, argv+argc); // passing pointers to argv values
    if(parseArgument(args, "--help") || parseArgument(args, "--h")){
        print_help();
        return 0;
    }
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
        std::cout << "Running without Instatniator, you can use program arguments: --s_type --s_max_it --phi \n\n";

        /*Read arguments*/
        SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
        int solver_max_it = parseArgument(args, "--s_max_it", 1000); //or 2e4;
        type_calc phi = parseArgument(args, "--phi", -600.0);


        
        /*Instantiate World*/ 
        world_ptr = std::make_unique<World>(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.2});
        world_ptr ->setTimeStart(); //optional, use it where you want your start of Wall time
        world_ptr->setTime(2e-11, 800);

        /*Instantiate objects*/ 
        world_ptr->addObject<Rectangle>(type_calc3(world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getX0()[2]),
         phi, type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], 0.05));
        
        world_ptr->addObject<Rectangle>(type_calc3(world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getXm()[2]),
         -phi, type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], 0.05));
        
        world_ptr->computeObjectID();


        type_calc E_ion_O = 1313.9*1000/Const::N_a; //in Joules
        /*Instantiate species*/ 
        species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 5e6, E_ion_O);
        species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 5e6);  // ions //last is mpw0
        species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e2);
        species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr, 5e6);      // electrons

        const type_calc num_den_neutrals = 1e17; 
        const type_calc num_den_ions = 1e10; //mean ion density
        const type_calc num_den_electrons = 1e15; //mean electron density
        type_calc T = 300;

        type_calc side_z = world_ptr->getXm()[2] - world_ptr->getX0()[2] - 0.1;
        species[0].loadParticleBoxThermal(world_ptr->getXc(), type_calc3(world_ptr->getL()[0]*0.3, world_ptr->getL()[1]*0.3, side_z), num_den_neutrals, T);
        species[3].loadParticleBoxThermal(type_calc3{world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getX0()[2] + 0.04}, type_calc3(world_ptr->getDx()[0], world_ptr->getDx()[1], world_ptr->getDx()[2]), num_den_electrons, T);
        
        /* Instantiate interactions*/
        type_calc rate_0_1 = 1e-2;
        type_calc rate_1_2 = 1e-3;
        interactions.emplace_back(std::make_unique<DSMC_MEX_IONIZATION>(species[0], species[1], species[3], *world_ptr));
        //interactions.emplace_back(std::make_unique<DSMC_MEX_IONIZATION>(species[1], species[2], species[3], *world_ptr));

        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[0], species[1], species[3], *world_ptr, rate_0_1));
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[1], species[2], species[3], *world_ptr, rate_1_2));

        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));


        /*Instantiate solver*/ 
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        solver_ptr->setReferenceValues(0, 1.5, num_den_ions);

        std::cout << "Solver Type: " << solver_type << "\n";
	    std::cout << "GS Solver max iterations: " << solver_ptr->get_GS_max_it() <<"\n";
        if(solver_type==PCG){
    	    std::cout << "PCG Solver max iterations: " << solver_ptr->get_PCG_max_it() <<"\n";
        }
	    std::cout << "Potential at one electrode: " << phi << " V" << std::endl;

    }
    Species& neutral_oxygen = species[0];
    Species& one_oxygen = species[1];
    Species& two_oxygen = species[2];
    Species& electrons = species[3];


    //////////////
    
    for (Species &sp : species) {
        std::cout << sp.name << " has " << sp.getNumParticles() << " particles" << std::endl;
    }

    //Output::fields(*world_ptr, species);
    solver_ptr->solve();
    solver_ptr->computeEF();
    Output::fields(*world_ptr, species);

    //world_ptr->setTimeStart();
    while(world_ptr->advanceTime()){

        for(std::unique_ptr<Source>& source: sources){
            source->sample();
        }

        for(std::unique_ptr<Interaction>& interaction: interactions){
            interaction->apply(world_ptr->getDt());
        }

        auto time_start = world_ptr->getWallTime();
        for(Species&  sp: species){
            sp.advance(neutral_oxygen, neutral_oxygen); //TODO correct it so for species neutrals refer to species neutrals
            sp.computeNumberDensity();
            sp.sampleMoments();
            sp.computeMacroParticlesCount();
        }

        world_ptr->computeChargeDensity(species);
        solver_ptr->solve();
        solver_ptr->computeEF();

        if(world_ptr->getTs()>5){
            for(Species&  sp: species){
                sp.updateAverages();
            }
        }
        
        Output::diagOutput(*world_ptr, species);
        Output::screenOutput(*world_ptr, species);
        int ts = world_ptr->getTs();
        if(ts%10 == 0 || world_ptr->isLastTimeStep()){ //|| (ts > 140 && ts < 160)){
            Output::fields(*world_ptr, species);
            //Output::particles(*world_ptr, species, 10000);
            std::cout << "Time taken so far: " << world_ptr->getWallTime() << std::endl;

        }
       
    }
    std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
    return 0;
}

