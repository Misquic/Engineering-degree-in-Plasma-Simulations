#include <iostream>
#include <memory>
#include "funkc.h"
#include "all.h"
#include "World.h"
#include "Field.h"
#include "Vec3.h"
#include "Outputs.h"
#include "PotentialSolver.h"
#include "Species.h"
#include "Object.h"
#include "Objects.h"

int main(int argc, char* argv[] ){

    std::vector<Species> species;
    //std::vector<std::unique_ptr<Object>> objects_ptrs; // using vector to object pointers so that we can use polimorphysm and store multiple different Objects in one container
    std::unique_ptr<Objects> objects_ptr;
    std::unique_ptr<PotentialSolver> solver_ptr;
    std::unique_ptr<World> world_ptr;

    if(argc > 1 && std::string(argv[1]) == "g"){
        //doesn't support Objects yet
        Instantiator instantiator(world_ptr, species, solver_ptr);
    }
    else{
        // Instantiate World
        world_ptr = std::make_unique<World>(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
        world_ptr->setTime(1e-7, 400);

        // Instantiate sphere
        type_calc phi_sphere = -100;
        if(argc > 1){
            phi_sphere = std::stof(argv[1]);
        }
        objects_ptr = std::make_unique<Objects>(*world_ptr);
        objects_ptr->addObject<Sphere>(type_calc3(0, 0, 0.15), phi_sphere, 0.05);

        // Instantiate species
        species.reserve(2);
        species.emplace_back("O+", 16 * Const::amu, Const::q_e, *world_ptr);  // ions
        species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr);      // electrons

        int3 num_ions_grid = {41, 41, 41};
        int3 num_eles_grid = {21, 21, 21};
        species[0].loadParticleBoxQS(world_ptr->getX0(), world_ptr->getXm(), 1e11, num_ions_grid);
        species[1].loadParticleBoxQS(world_ptr->getX0(), world_ptr->getXc(), 1e11, num_eles_grid);
    }

    //////////////
    
    for (Species &sp : species) {
        std::cout << sp.name << " has " << sp.getNumParticles() << " particles" << std::endl;
    }

    // Instantiate solver
    if (!solver_ptr) {
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, 1e4, 1e-4);  // Jeśli solver nie został zainicjalizowany
    }

    solver_ptr->solveGS();
    solver_ptr->computeEF();

    world_ptr->setTimeStart();
    while(world_ptr->advanceTime()){

        for(Species&  sp: species){
            sp.advance();
            sp.computeNumberDensity();
        }

        world_ptr->computeChargeDensity(species);
        solver_ptr->solveGS();
        solver_ptr->computeEF();

        Output::diagOutput(*world_ptr, species);

        if(world_ptr->getTs()%100 == 0 || world_ptr->isLastTimeStep()){
            Output::screenOutput(*world_ptr, species);
            Output::fields(*world_ptr, species);
        }
       
    }
    std::cout << "Simulation took " << world_ptr->getWallTime() << "seconds" << std::endl;
    return 0;
}

