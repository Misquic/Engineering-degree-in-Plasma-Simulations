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

//tags: flow aroung a shape(sphere), shape inside, inlet, neuman elsewhere, Maxwellian velocity distribution//

int main(int argc, char* argv[] ){

    std::vector<Species> species;
    std::vector<ColdBeamSource> sources;
    //std::vector<std::unique_ptr<Object>> objects_ptrs; // using vector to object pointers so that we can use polimorphysm and store multiple different Objects in one container
    std::unique_ptr<PotentialSolver> solver_ptr;
    std::unique_ptr<World> world_ptr;

    if(argc > 1 && std::string(argv[1]) == "g"){
        //doesn't support Objects yet
        Instantiator instantiator(world_ptr, species, solver_ptr);
        //set ref values solver
        //ad objects and calc object_id
    }
    else{
        std::cout << "else\n";
        // Instantiate World
        world_ptr = std::make_unique<World>(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
        world_ptr->setTime(1e-7, 400);

        // Instantiate sphere
        type_calc phi_sphere = -100;
        if(argc > 1){
            phi_sphere = std::stof(argv[1]);
        }
	    std::cout<<"Sphere potential: "<<phi_sphere<<" V"<<std::endl;

        world_ptr->addObject<Sphere>(type_calc3(0, 0, 0.15), phi_sphere, 0.05);
        //world_ptr->addObject<Rectangle>(type_calc3(0, 0, 0.15), phi_sphere, type_calc3(0.1, 0.07, 0.2), type_calc3(1, 1, 1));
        world_ptr->computeObjectID();
        std::string inlet_Face = "z-";
        world_ptr->addInlet(inlet_Face);

        // Instantiate species
        species.reserve(2);
        species.emplace_back("O+", 16 * Const::amu, Const::q_e, *world_ptr, 5e1);  // ions //last is mpw0
        //species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr);      // electrons

        // Instantiate sources
        const type_calc ndi = 1e10; //mean ion density
        sources.reserve(1);
        sources.emplace_back(species[0], *world_ptr, 7000, ndi);  // &sp, &world, v_drift, den 
        

        // Instantiate solver
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr,20000,1e-4,SolverType::PCG);  // Jeśli solver nie został zainicjalizowany
        solver_ptr->setReferenceValues(0, 1.5, ndi);
    }

    //////////////
    
    for (Species &sp : species) {
        std::cout << sp.name << " has " << sp.getNumParticles() << " particles" << std::endl;
    }

    // Instantiate solver

    //Output::fields(*world_ptr, species);
    solver_ptr->solve();
    solver_ptr->computeEF();
    Output::fields(*world_ptr, species);

    world_ptr->setTimeStart();
    while(world_ptr->advanceTime()){

        for(ColdBeamSource& source: sources){
            source.sample();
        }
        for(Species&  sp: species){
            sp.advance();
            sp.computeNumberDensity();
        }

        world_ptr->computeChargeDensity(species);
        solver_ptr->solve();
        solver_ptr->computeEF();

        if(world_ptr->steadyState(species)){
            for(Species&  sp: species){
                sp.updateAverages();
            }
        }

        //Output::diagOutput(*world_ptr, species);
        Output::screenOutput(*world_ptr, species);

        if(world_ptr->getTs()%20 == 0 || world_ptr->isLastTimeStep()){
            Output::fields(*world_ptr, species);
            std::cout << "Time taken so far: " << world_ptr->getWallTime() << std::endl;
        }
       
    }
    std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
    return 0;
}

