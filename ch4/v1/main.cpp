// #include <iostream>
// #include <memory>
// #include "funkc.h"
// #include "Instantiator.h"
// #include "all.h"
// #include "World.h"
// #include "Field.h"
// #include "Vec3.h"
// #include "Outputs.h"
// #include "PotentialSolver.h"
// #include "Species.h"
// #include "Object.h"
// #include "Source.h"
// #include "Rnd.h"
// #include "Interactions.h"

// /*tags: flow aroung a shape(sphere and quboid), shape inside, inlet, neuman elsewhere, Maxwellian velocity distribution*/
// int main(int argc, char* argv[] ){
//     // TESTS: ////////////////////////////


//     // World world(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
//     // Species sp("Sp", Const::amu, Const::q_e, world, 200);

//     // type_calc T = 300, vel = 1e3, mpw = 200;
//     // sp.sampleVthVariableMpw(T, vel, mpw);
//     // return 0;

//     ///////////////////////////////////////

//     std::vector<std::string> args(argv + 1, argv+argc); // passing pointers to argv values
//     if(parseArgument(args, "--help") || parseArgument(args, "--h")){
//         print_help();
//         return 0;
//     }
//     std::vector<Species> species; //species container
//     std::vector<std::unique_ptr<Source>> sources; //sources container
//     std::unique_ptr<PotentialSolver> solver_ptr;
//     std::unique_ptr<World> world_ptr;
//     std::vector<std::unique_ptr<Interaction>> interactions;

//     if(parseArgument(args, "--i")){
//         /*doesn't support Objects yet
//           TODO trash it? command line arguments are usable and better, but Instatniator might be better for future UI?*/
//         Instantiator instantiator(world_ptr, species, solver_ptr);

//     }
//     else{
//         std::cout << "Running without Instatniator, you can use program arguments: --s_type --s_max_it --sphere_phi \n\n";

//         /*Read arguments*/
//         SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
//         int solver_max_it = parseArgument(args, "--s_max_it", 1000); //or 2e4;
//         type_calc phi_sphere = parseArgument(args, "--sphere_phi", -100.0);


        
//         /*Instantiate World*/ 
//         // world_ptr = std::make_unique<World>(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
//         world_ptr = std::make_unique<World>(41, 21, 41, type_calc3{-0.2,-0.1,0}, type_calc3{0.2,0.1,0.4});
//         world_ptr ->setTimeStart();
//         world_ptr->setTime(1e-7, 600);

//         /*Instantiate sphere*/ 
//         world_ptr->addObject<Sphere>(type_calc3(0, 0, 0.15), phi_sphere, 0.05);
//         // world_ptr->addObject<Rectangle>(type_calc3(0, 0, 0.35), phi_sphere, type_calc3(0.05, 0.05, 0.05));
//         world_ptr->computeObjectID();
//         std::string inlet_Face = "-z";
//         world_ptr->addInlet(inlet_Face);

//         /*Instantiate species*/ 
//         // species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 1e2);  // ions //last is mpw0
//         // species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e1);
//         // species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 2e3);
//         species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 1e12);
//         /*spherium (for rectangle also name temporary)*/
//         // species.emplace_back("Sph", 100*Const::amu, 0, *world_ptr, 2e2);
//         // species.emplace_back("Sph+", 100*Const::amu, Const::q_e, *world_ptr, 2e2);
//         // species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr);      // electrons


//         /* Instantiate sources */
//         const type_calc num_den_neutrals = 2e20; 
//         const type_calc num_den_ions = 1e10; //mean ion density
//         type_calc T = 0;
//         // sources.emplace_back(species[0], *world_ptr, 7000, 0.8*num_den_ions, T, inlet_Face); //O+
//         // sources.emplace_back(species[1], *world_ptr, 7000, 0.1*num_den_ions, T, inlet_Face); //O++
//         // sources.emplace_back(std::make_unique<ColdBeamSource>(species[0], *world_ptr, 7000, 0.8*num_den_ions, inlet_Face)); //O+
//         // sources.emplace_back(std::make_unique<ColdBeamSource>(species[1], *world_ptr, 7000, 0.1*num_den_ions, inlet_Face)); //O++
//         sources.emplace_back(std::make_unique<WarmBeamSource>(species[0], *world_ptr, 7000, num_den_neutrals, T, inlet_Face)); //O
//         //sources.emplace_back(std::make_unique<WarmBeamSource>(species[1], *world_ptr, 7000, 0.1*num_den_ions, T, inlet_Face)); //O++
//         //sources.emplace_back(species[2], *world_ptr, 7000, num_den_neutrals); //O
        
//         /* Instantiate interactions*/
//         // type_calc rate = 1e-4;
//         // interactions.emplace_back(std::make_unique<ChemistryIonize>(species[3], species[4], *world_ptr, rate));
//         interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));


//         /*Instantiate solver*/ 
//         solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
//         solver_ptr->setReferenceValues(0, 1.5, num_den_ions);

//         std::cout << "Solver Type: " << solver_type << "\n";
// 	    std::cout << "GS Solver max iterations: " << solver_ptr->get_GS_max_it() <<"\n";
//         if(solver_type==PCG){
//     	    std::cout << "PCG Solver max iterations: " << solver_ptr->get_PCG_max_it() <<"\n";
//         }
// 	    std::cout << "Sphere potential: " << phi_sphere << " V" << std::endl;

//     }

//     // Species& neutrals = species[2];
//     // Species& spherium = species[3];
//     Species& neutrals = species[0];


//     //////////////
    
//     for (Species &sp : species) {
//         std::cout << sp.name << " has " << sp.getNumParticles() << " particles" << std::endl;
//     }

//     //Output::fields(*world_ptr, species);
//     solver_ptr->solve();
//     solver_ptr->computeEF();
//     Output::fields(*world_ptr, species);

//     //world_ptr->setTimeStart();
//     while(world_ptr->advanceTime()){

//         for(std::unique_ptr<Source>& source: sources){
//             source->sample();
//         }

//         for(std::unique_ptr<Interaction>& interaction: interactions){
//             interaction->apply(world_ptr->getDt());
//         }

//         auto time_start = world_ptr->getWallTime();
//         for(Species&  sp: species){
//             sp.advance(neutrals, neutrals); //TODO correct it so for species neutrals refer to species neutrals
//             sp.computeNumberDensity();
//             sp.sampleMoments();
//             sp.computeMacroParticlesCount();
//         }

//         world_ptr->computeChargeDensity(species);
//         //solver_ptr->solve();
//         //solver_ptr->computeEF();

//         if(world_ptr->getTs()>5){
//             for(Species&  sp: species){
//                 sp.updateAverages();
//             }
//         }
        
//         Output::diagOutput(*world_ptr, species);
//         Output::screenOutput(*world_ptr, species);
//         int ts = world_ptr->getTs();
//         if(ts%10 == 0 || world_ptr->isLastTimeStep()){ //|| (ts > 140 && ts < 160)){
//             Output::fields(*world_ptr, species);
//             //Output::particles(*world_ptr, species, 10000);
//             std::cout << "Time taken so far: " << world_ptr->getWallTime() << std::endl;

//         }
       
//     }
//     std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
//     return 0;
// }

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

/*tags: flow aroung a shape(sphere and quboid), shape inside, inlet, neuman elsewhere, Maxwellian velocity distribution*/
int main(int argc, char* argv[] ){
    // TESTS: ////////////////////////////


    // World world(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
    // Species sp("Sp", Const::amu, Const::q_e, world, 200);

    // type_calc T = 300, vel = 1e3, mpw = 200;
    // sp.sampleVthVariableMpw(T, vel, mpw);
    // return 0;

    ///////////////////////////////////////

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
        std::cout << "Running without Instatniator, you can use program arguments: --s_type --s_max_it --sphere_phi \n\n";

        /*Read arguments*/
        SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
        int solver_max_it = parseArgument(args, "--s_max_it", 1000); //or 2e4;
        type_calc phi_sphere = parseArgument(args, "--sphere_phi", -100.0);


        
        /*Instantiate World*/ 
        world_ptr = std::make_unique<World>(21, 21, 41, type_calc3{-0.1, -0.1, 0.0}, type_calc3{0.1, 0.1, 0.4});
        // world_ptr = std::make_unique<World>(41, 21, 41, type_calc3{-0.2,-0.1,0}, type_calc3{0.2,0.1,0.4});
        world_ptr ->setTimeStart();
        world_ptr->setTime(1e-7, 600);

        /*Instantiate sphere*/ 
        world_ptr->addObject<Sphere>(type_calc3(0, 0, 0.15), phi_sphere, 0.05);
        // world_ptr->addObject<Rectangle>(type_calc3(0, 0, 0.35), phi_sphere, type_calc3(0.05, 0.05, 0.05));
        world_ptr->computeObjectID();
        std::string inlet_Face = "-z";
        world_ptr->addInlet(inlet_Face);

        /*Instantiate species*/ 
        species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 1e2);  // ions //last is mpw0
        species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e1);
        species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 2e3);
        // species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 1e12);
        /*spherium (for rectangle also name temporary)*/
        species.emplace_back("Sph", 100*Const::amu, 0, *world_ptr, 2e2);
        species.emplace_back("Sph+", 100*Const::amu, Const::q_e, *world_ptr, 2e2);
        // species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr);      // electrons


        /* Instantiate sources */
        // const type_calc num_den_neutrals = 2e20; 
        // const type_calc num_den_ions = 1e10; //mean ion density
        const type_calc num_den_neutrals = 1e13; 
        const type_calc num_den_ions = 1e10; //mean ion density
        type_calc T = 300;
        // sources.emplace_back(species[0], *world_ptr, 7000, 0.8*num_den_ions, T, inlet_Face); //O+
        // // sources.emplace_back(species[1], *world_ptr, 7000, 0.1*num_den_ions, T, inlet_Face); //O++
        // sources.emplace_back(std::make_unique<ColdBeamSource>(species[0], *world_ptr, 7000, 0.8*num_den_ions, inlet_Face)); //O+
        // sources.emplace_back(std::make_unique<ColdBeamSource>(species[1], *world_ptr, 7000, 0.1*num_den_ions, inlet_Face)); //O++
        sources.emplace_back(std::make_unique<WarmBeamSource>(species[0], *world_ptr, 7000, 0.8*num_den_ions, T, inlet_Face)); //O+
        sources.emplace_back(std::make_unique<WarmBeamSource>(species[1], *world_ptr, 7000, 0.1*num_den_ions, T, inlet_Face)); //O++
        //sources.emplace_back(species[2], *world_ptr, 7000, num_den_neutrals); //O
        
        /* Instantiate interactions*/
        type_calc rate = 1e-4;
        interactions.emplace_back(std::make_unique<ChemistryIonize>(species[3], species[4], *world_ptr, rate));
        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));


        /*Instantiate solver*/ 
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        solver_ptr->setReferenceValues(0, 1.5, num_den_ions);

        std::cout << "Solver Type: " << solver_type << "\n";
	    std::cout << "GS Solver max iterations: " << solver_ptr->get_GS_max_it() <<"\n";
        if(solver_type==PCG){
    	    std::cout << "PCG Solver max iterations: " << solver_ptr->get_PCG_max_it() <<"\n";
        }
	    std::cout << "Sphere potential: " << phi_sphere << " V" << std::endl;

    }

    // Species& neutrals = species[2];
    // Species& spherium = species[3];
    Species& neutrals = species[2];
    Species& spherium = species[3];


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
            sp.advance(neutrals, spherium); //TODO correct it so for species neutrals refer to species neutrals
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

