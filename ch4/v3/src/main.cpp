#include <iostream>
#include <memory>
#include <thread>
#include "funkc.h"
// #include "Instantiator.h"
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
#include "Config.h"

/*TODO: add option to turn off spherium, improve ionisation, add space varying epsilon, add instantiator back, add loading with velocity sampling*/
int main(int argc, char* argv[] ){

    // Config::getInstance().setNUM_THREADS(10);
    // std::cout << splitIntoChunks(1)<< "\n";
    // std::cout << splitIntoChunks(9)<< "\n";
    // std::cout << splitIntoChunks(10)<< "\n";
    // std::cout << splitIntoChunks(11)<< "\n";
    // std::cout << splitIntoChunks(12)<< "\n";

    // return 0;

    std::ofstream errorFile("outputs/output.txt"); // File to writing errors
    if (!errorFile) {
        std::cerr << "Cannot open file to write std::cerr !" << std::endl;
        return 1;
    }

    std::cerr.rdbuf(errorFile.rdbuf()); 

    #ifdef DEBUG
    std::cout << "DEBUG mode\n";
    dmsg("Debug msg\n");
    #else
    std::cout << "RELEASE mode\n";
    #endif

    std::vector<std::string> args(argv + 1, argv+argc); // passing pointers to argv values
    if(parseArgument(args, "--help") || parseArgument(args, "--h")){
        print_help();
        return 0;
    }
    
    Config& config = Config::getInstance();
    config.setSUBCYCLING(parseArgument(args, "--subcycling", true));
    config.setMULTITHREADING(parseArgument(args, "--multithreading", true));
    config.setNUM_THREADS(parseArgument(args, "--num_threads", (unsigned int)15));
    config.setMERGING(parseArgument(args, "--merging", true));

    Output::modes OUTPUT = parseArgument(args, "--output", Output::modes::fields); // 0 none // 1 all // 2 only screen on // 3 only screen and fields on // 4 only screen, fields and particles on etc.

    std::vector<Species> species; //species container
    std::vector<std::unique_ptr<Source>> sources; //sources container
    std::unique_ptr<PotentialSolver> solver_ptr;
    std::unique_ptr<World> world_ptr;
    std::vector<std::unique_ptr<Interaction>> interactions; //interactions container

    {
        std::cout << "Running without Instatniator, you can use program arguments: --subcycling [<bool>] --s_type [PCG, GS] --s_max_it [<int>] --phi [<double>] --num_ts [<int>] --dt [<double>] \n\n";

        /*Read arguments*/
        SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
        int solver_max_it = parseArgument(args, "--s_max_it", 8000); //or 2e4;
        type_calc solver_tol = parseArgument(args, "--s_tol", 1); //or 1e-4;
        type_calc phi = parseArgument(args, "--phi", -4000.0);
        int num_ts = parseArgument(args, "--num_ts", (int)4000);
        type_calc dt = parseArgument(args, "--dt", 10e-13);
        
        /*Instantiate World*/ 
        world_ptr = std::make_unique<World>(41, 41, 61, type_calc3{-0.004, -0.004, 0.0}, type_calc3{0.004, 0.004, 0.005}); //0.4cm x 0.4cm x 0.5 cm
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
        species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 5e11, E_ion_O);
        species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 100);  // ions //last is mpw0
        // species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e0);
        species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr, 100);      // electrons

        const type_calc num_den_neutrals = 1e25; // final should be 1e25 
        const type_calc num_den_ions = 1e10; //mean ion density
        const type_calc num_den_electrons = 1e16; //mean electron density
        type_calc T = 300;
        type_calc T_eles = 10*T;

        // type_calc side_z = (world_ptr->getXm()[2] - world_ptr->getX0()[2])*0.9;
        species[2].loadParticleBoxThermal(type_calc3{world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getXc()[2] - 0.5*world_ptr->getL()[2]*0.85}, type_calc3(world_ptr->getL()[0]*0.01, world_ptr->getL()[1]*0.01, world_ptr->getL()[2]*0.02), num_den_electrons, T_eles);

        species[0].loadParticleBoxThermal(world_ptr->getXc(), type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], world_ptr->getL()[2]*0.9), num_den_neutrals, T);
        
        /* Instantiate interactions*/
        interactions.emplace_back(std::make_unique<MC_MEX_Ionization>(species[0], species[1], species[2], *world_ptr));
        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[2], *world_ptr));
        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[1], *world_ptr));

        // interactions.emplace_back(std::make_unique<DSMC_MEX_Ionization>(species[1], species[2], species[3], *world_ptr));

        // type_calc rate_0_1 = 1e-2;
        // type_calc rate_1_2 = 1e-3;
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[0], species[1], species[3], *world_ptr, rate_0_1));
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[1], species[2], species[3], *world_ptr, rate_1_2));

        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));

        /*Instantiate solver*/ 
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it, solver_tol, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        // solver_ptr->setReferenceValues(0, num_den_ions, T); //simulationg electrons via boltzman relationship
        solver_ptr->setReferenceValues(0, 0, 1e20); //simulationg electrons directly 0 in ne zeroes infuence from boltzman relationship


        std::cout << "\nOptions:\n";
        std::cout << "Sybcycling: " << config.getSUBCYCLING() << "\n";
        std::cout << "Multithreading: " << config.getMULTITHREADING() << "\n";
        if(config.getMULTITHREADING()){
            std::cout << "Num Threads: " << config.getNUM_THREADS() << "\n";
        }
        std::cout << "Merging " << config.getMERGING() << "\n";
        std::cout << "output mode: " << OUTPUT << "\n";
        std::cout << "Solver Type: " << solver_type << "\n";
	    std::cout << "GS Solver max iterations: " << solver_ptr->get_GS_max_it() <<"\n";
        if(solver_type==PCG){
    	    std::cout << "PCG Solver max iterations: " << solver_ptr->get_PCG_max_it() <<"\n";
        }
	    std::cout << "Potential at one electrode: " << phi << " V\n" << std::endl;
	    std::cout << "Time step: " << dt << " s\n" << std::endl;
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

    // Output::fields(*world_ptr, species);
    solver_ptr->solve();
    solver_ptr->computeEF();
    Output::fieldsOutput(*world_ptr, species);

    // world_ptr->setTimeStart();
    type_calc dt = world_ptr->getDt();
    while(world_ptr->advanceTime()){
        dmsg("\ntime step: " << world_ptr->getTs() << "\n");
        if(world_ptr->getTs() > 250 && world_ptr->getTs()%50 == 0 && config.getMERGING()){
            if(!config.getMULTITHREADING()){
                for(Species&  sp: species){ 
                    sp.merge();
                }
            }
            else{
                std::vector<std::thread> threads_merge;
                for(int i = 0; i < species.size(); i++){
                    threads_merge.emplace_back(std::thread(&Species::merge, &species[i]));
                    // sp.merge();
                }
                for(auto &th: threads_merge){
                    th.join();
                }
            }
        }
        dmsg("source start\n");
        for(std::unique_ptr<Source>& source: sources){
            source->sample();
        }
        dmsg("interaction start\n");
        for(std::unique_ptr<Interaction>& interaction: interactions){
            interaction->apply(world_ptr->getDt());
        }

        dmsg("subcycling start\n");
        ///////////////////////////////////////
        // for(Species& sp: species){ // make it look like this!
        //     sp.advance();
        // }
        ///////////////////////////////////////
        if(config.getSUBCYCLING() == true){
            for(Species&  sp: species){// TODO correct it so for species neutrals refer to species neutrals
                if(sp.name == electrons.name){ // implemented subcycling //electrons
                    dmsg("electrons start\n");
                    type_calc time_start = world_ptr->getWallTime();
                    sp.advanceElectrons(dt);
                    std::cout << "\ntime electrons [s] : " << world_ptr->getWallTime()-time_start;

                    sp.sampleMoments();
                    sp.computeNumberDensity();
                    sp.computeMacroParticlesCount();
                
                }
                else if(sp.charge!= 0 && world_ptr->getTs()%10 == 0){ //ions
                    dmsg("ions start\n");
                    sp.advanceNonElectron(neutral_oxygen, neutral_oxygen, dt*10); //sputtering according to config
                    sp.computeNumberDensity();
                    sp.sampleMoments();
                    sp.computeMacroParticlesCount();
                }
                else if(world_ptr->getTs()%100 == 0){ // neutrals 
                    sp.advanceNonElectron(neutral_oxygen, neutral_oxygen, dt*50); //sputtering according to config, neutrals can't do it
                    sp.computeNumberDensity();
                    sp.sampleMoments();
                    sp.computeMacroParticlesCount();
                }
                // sp.advance(neutral_oxygen, neutral_oxygen, dt); //normal
            }
        }
        else{
            for(Species&  sp: species){// TODO correct it so for species neutrals refer to species neutrals
                if(sp.name == electrons.name){ // implemented subcycling //electrons
                    sp.advanceElectrons(dt);
                }
                else{ //ions
                    sp.advanceNonElectron(neutral_oxygen, neutral_oxygen, dt); // sputtering according to config
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

        dmsg("density start\n");
        world_ptr->computeChargeDensity(species);
        // solver_ptr->solve();
        // solver_ptr->computeEF();


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
        //std::cout << std::flush <<"\n";

       
    }
    if(world_ptr->getWallTime()>60){
        std::cout << "Simulation took " << world_ptr->getWallTime()/60 << " minutes." << std::endl;
    }else{
        std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
    }
    return 0;
}
