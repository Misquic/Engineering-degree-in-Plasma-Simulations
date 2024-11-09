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

//tags: flow aroung a shape(sphere), shape inside, inlet, neuman elsewhere, Maxwellian velocity distribution//
int main(int argc, char* argv[] ){

    // World world(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
    // world.addObject<Rectangle>(type_calc3(0, 0, 0), -100,type_calc3{0.1, 0.2, 0.3});
    // // world.addObject<Sphere>(type_calc3(0, 0, 0.15), -100, 0.05);

    // type_calc3 x1{0.2, -0.39, -0.6}, x2{0.0, 0.0, 0.0};
    // type_calc tp{};
    // type_calc3 n{};
    // type_calc3 x3{};
    // if(world.inObject(x2) &&  !world.inObject(x1)){
    //     world.lineIntersect(x1, x2, 1, tp, x3, n);
    // }
    

    // std::cout << "tp: " << tp << "\n";
    // std::cout << std::setw(0) << x1 << " " << x2 << " "  << x3 << " " << n << "\n";

    // return 0;


    std::vector<std::string> args(argv + 1, argv+argc); // passing pointers to argv values
    if(parseArgument(args, "--help") || parseArgument(args, "--h")){
        print_help();
        return 0;
    }
    std::vector<Species> species;
    std::vector<ColdBeamSource> sources;
    //std::vector<std::unique_ptr<Object>> objects_ptrs; // using vector to object pointers so that we can use polimorphysm and store multiple different Objects in one container
    std::unique_ptr<PotentialSolver> solver_ptr;
    std::unique_ptr<World> world_ptr;

    if(parseArgument(args, "--i")){
        //doesn't support Objects yet
        //TODO trash it? command line arguments are usable and better, but Instatniator might be better for future UI?
        Instantiator instantiator(world_ptr, species, solver_ptr);
        //set ref values solver
        //ad objects and calc object_id
    }
    else{
        std::cout << "Running without Instatniator, you can use program arguments: --s_type --s_max_it --sphere_phi \n\n";

        //read arguments
        SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
        int solver_max_it = parseArgument(args, "--s_max_it", 1000); //or 2e4;
        type_calc phi_sphere = parseArgument(args, "--sphere_phi", -100.0);


        
        // Instantiate World
        world_ptr = std::make_unique<World>(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
        world_ptr ->setTimeStart();
        world_ptr->setTime(1e-7, 600);

        // Instantiate sphere
        world_ptr->addObject<Sphere>(type_calc3(0, 0, 0.15), phi_sphere, 0.05);
        //world_ptr->addObject<Rectangle>(type_calc3(0, 0, 0.15), phi_sphere, type_calc3(0.1, 0.07, 0.2), type_calc3(1, 1, 1));
        world_ptr->computeObjectID();
        std::string inlet_Face = "z-";
        world_ptr->addInlet(inlet_Face);

        // Instantiate species
        //species.reserve(3);
        species.emplace_back("O+", 16 * Const::amu, Const::q_e, *world_ptr, 1e2);  // ions //last is mpw0
        species.emplace_back("O++", 16 * Const::amu, 2*Const::q_e, *world_ptr, 5e1);
        species.emplace_back("O", 16 * Const::amu, 0, *world_ptr, 1e5);
        // species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr);      // electrons

        // Instantiate sources
        const type_calc num_den_ions = 1e10; //mean ion density
        const type_calc num_den_neutrals = 1e13; 
        //sources.reserve(3);
        sources.emplace_back(species[0], *world_ptr, 7000, 0.8*num_den_ions); //O+
        sources.emplace_back(species[1], *world_ptr, 7000, 0.1*num_den_ions); //O++
        sources.emplace_back(species[2], *world_ptr, 7000, num_den_neutrals); //O
        

        // Instantiate solver
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        solver_ptr->setReferenceValues(0, 1.5, num_den_ions);

        std::cout << "Solver Type: " << solver_type << "\n";
	    std::cout << "GS Solver max iterations: " << solver_ptr->get_GS_max_it() <<"\n";
        if(solver_type==PCG){
    	    std::cout << "PCG Solver max iterations: " << solver_ptr->get_PCG_max_it() <<"\n";
        }
	    std::cout << "Sphere potential: " << phi_sphere << " V" << std::endl;

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

    //world_ptr->setTimeStart();
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

        Output::diagOutput(*world_ptr, species);
        Output::screenOutput(*world_ptr, species);
        int ts = world_ptr->getTs();
        if(ts%10 == 0 || world_ptr->isLastTimeStep()){ //|| (ts > 140 && ts < 160)){
            Output::fields(*world_ptr, species);
            std::cout << "Time taken so far: " << world_ptr->getWallTime() << std::endl;
        }
       
    }
    std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
    return 0;
}

