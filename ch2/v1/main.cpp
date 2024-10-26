#include <iostream>
#include <chrono>
//#include <funkc.h>
#include "all.h"
#include "World.h"
#include "Field.h"
#include "Vec3.h"
#include "Outputs.h"
#include "PotentialSolver.h"
#include "Species.h"


//TODO dodać prosty interfejs żeby przyjemniej zmieniać dane, żeby nie kompilować ciągle
int main(){

    //instantiate world
    World world(21, 21, 21, {-0.1, -0.1, 0.0}, {0.1, 0.1, 0.2});
    world.setTime(2e-10, 1e4);

    //instantiate species
    std::vector<Species> species;
    species.reserve(2);
    species.emplace_back("O+", 16*Const::amu, Const::q_e, world);  //ions
    species.emplace_back("e-", Const::m_e, -Const::q_e, world);    //electrons



    int3 num_ions_grid = {41, 41, 41};
    int3 num_eles_grid = {21, 21, 21};
    species[0].loadParticleBoxQS(world.getX0(), world.getXm(), 1e11, num_ions_grid);
    species[1].loadParticleBoxQS(world.getX0(), world.getXc(), 1e11, num_eles_grid);
    // species[0].loadParticleBox(world.getX0(), world.getXm(), 1e11, 41*41*41);
    // species[1].loadParticleBox(world.getX0(), world.getXc(), 1e11, 21*21*21);
    
    for ( Species &sp : species ){
        std::cout<< sp.name << " has " << sp.getNumParticles() << " particles" << std::endl;
    }

    //instantiate solver
    PotentialSolver solver(world, 1e4, 1e-4);
    solver.solveGS();
    solver.computeEF();

    world.setTimeStart();
    while(world.advanceTime()){
///////////////////////////////////////


        species[0].advance();
        species[1].advance();
        //Output::fields(world, species, "advance");

        species[0].computeNumberDensity();
        species[1].computeNumberDensity();
        //Output::fields(world, species, "numdens");

        world.computeChargeDensity(species);
        //Output::fields(world, species, "chargedens");

        solver.solveGS();
        //Output::fields(world, species, "GS");

        solver.computeEF();
        //Output::fields(world, species, "EF");

        Output::diagOutput(world, species);

        //Output::fields(world, species);


        // Output::fields(world, species);//

        if(world.getTs()%100 == 0 || world.isLastTimeStep()){
            Output::screenOutput(world, species);
            Output::fields(world, species);
        }
       
    }
    std::cout << "Simulation took " << world.getWallTime() << "seconds" << std::endl;
    return 0;
}

