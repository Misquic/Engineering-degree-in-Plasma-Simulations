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

// #define OSCILATIONS
// #define DISCHARGE
#define SPHERE

#ifdef OSCILATIONS
int main(int argc, char* argv[]){

    std::vector<std::string> args(argv + 1, argv+argc); // passing pointers to argv values
    type_calc dx = parseArgument(args, "--dx", 0.01);
    bool QS = parseArgument(args, "--QS", false);
    bool T = parseArgument(args, "--T", 0);
    type_calc dt = parseArgument(args, "--dt", 2e-10);
    int num_ts = parseArgument(args, "--num_ts", 1e5);
    int solver_max_it = parseArgument(args, "--s_max_it", 2000); //or 2e4;
    SolverType solver_type = parseArgument(args, "--s_type", SolverType::PCG);
    int W = parseArgument(args, "--W", 1);
    type_calc den = parseArgument(args, "--den", 1e11);
    bool SUBCYCLING = parseArgument(args, "--subcycling", false);
    Output::modes OUTPUT = parseArgument(args, "--output", Output::modes::all); // 0 none // 1 all // 2 only screen on // 3 only screen and fields on // 4 only screen, fields and particles on etc.

    std::vector<Species> species; //species container
    std::vector<std::unique_ptr<Source>> sources; //sources container
    std::vector<std::unique_ptr<Interaction>> interactions;

    int3 grid(61, 31, 31);
    World world(grid);
    world.setExtents({0.0, 0.0, 0.0}, {0.2, 0.1, 0.1});
    world.dirichletBoundries();
    world.setTime(dt, num_ts);
    type_calc lambda = parseArgument(args, "--lambda", world.getL()[0]/4);
    type_calc A = parseArgument(args, "--A", 0.005);

    species.emplace_back("O+", 16*Const::amu, Const::q_e, world, 5*W);  //ions
    species.emplace_back("e-", Const::m_e, -Const::q_e, world, W);    //electrons
    species.emplace_back("e-_sample", Const::m_e, -Const::q_e, world, 1);    //electrons
    Rectangle box(world.getXc(), 0, world.getL()*0.2);

    int3 num_ions_grid = {101, 101, 101};
    int3 num_eles_grid = {101, 101, 101};
    
    if(QS){
        species[1].loadParticleSphereThermal(world.getXc()+type_calc3(world.getXc()[0]*dx,0,0), world.getL()[0]/4, den, T);
        species[0].loadParticleBoxQS(world.getXc(), world.getL(), den, num_ions_grid);
    }
    else{
        // // species[0].loadParticleBoxThermal(world.getXc(), world.getL(), den, T);
        // // species[1].loadParticleBoxThermal(world.getXc()+type_calc3(world.getXc()[0]*dx,0,0), world.getL()/2, den, T);        
        // // species[1].loadParticleSphereThermal(world.getXc()+type_calc3(world.getXc()[0]*dx,0,0), (world.getL()[0]/4), den, T);   
        // species[0].loadParticleBoxThermal(world.getXc(), world.getL(), den, T);
        // std::cout << "wtf\n";

        // species[1].loadLangumir(lambda, A, den, T);      
        species[0].loadParticleBoxThermal(world.getXc(), world.getL().elWiseMult({0.90, 1,1}), den, T);
        species[1].loadParticleBoxThermal(world.getXc(), world.getL().elWiseMult({0.90, 1,1}), den, T);
        species[2].loadParticleBoxThermal(world.getXc(), world.getL().elWiseMult({0.90, 1,1}), 1000/world.getCellVolume()/world.getNumCells(), T);
        // species[1].move(lambda, A);
        species[1].move({dx,0,0});
        species[2].move({dx,0,0});


        // species[0].loadParticleBoxThermal(world.getXc(), world.getL()*0.9, den, T);
        // species[1].loadParticleBoxThermal(world.getXc(), world.getL()*0.9, den, T);

        std::cout << "częstość langmuira: " << std::sqrt((Const::q_e*Const::q_e*den)/(Const::eps_0*Const::m_e)) << "\n";
    }

    PotentialSolver solver(world, solver_max_it, 1e-4, solver_type);
    solver.setReferenceValues(0, 0, 1e20); //simulationg electrons directly 0 in ne zeroes infuence from boltzman relationship

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for ( Species &sp : species ){
        std::cout<< sp.name << " has " << sp.getNumParticles() << " particles with weight " << sp.getPartRef()[0].macro_weight << std::endl;
    }

    world.setTimeStart(); //optional, use it where you want your start of Wall time

    solver.solveGS();
    solver.computeEF();
    Output::fieldsOutput(world, species);
    bool first = true;
    while(world.advanceTime()){

        // if(world.steadyState() && first == true){
        //     first = false;
        //     species[1].move(lambda, A);
        //     std::cout << "moving particles: " << type_calc3{A,0,0}.elWiseMult(world.getL());
        // }

        for(std::unique_ptr<Source>& source: sources){
            source->sample();
        }

        for(std::unique_ptr<Interaction>& interaction: interactions){
            interaction->apply(world.getDt());
        }

        for(Species&  sp: species){
            if(sp.charge < 0){
                sp.advanceElectrons(dt);
            }
            else if(world.getTs()%100){
                sp.advanceElectrons(dt*100);
            }

            sp.sampleMoments();
            sp.computeNumberDensity();
            sp.computeMacroParticlesCount();
        }


        if(world.getTs()>5){ //moved before next section
            for(Species&  sp: species){
                sp.updateAverages();
            }
        }

        world.computeChargeDensity(species);
        solver.solve();
        solver.computeEF();


        if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::diagnostics ){
            Output::diagOutput(world, species);
        }

        int ts = world.getTs();
        if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::screen ){
            Output::screenOutput(world, species);
            // std::cout << "\n";
        }
        if(ts%30 == 0 || world.isLastTimeStep()){ //|| (ts > 140 && ts < 160)){
            if(!world.steadyState() && ts>5){
                world.steadyState(species);
            }
            
            if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::fields ){
                Output::fieldsOutput(world, species);
            }
            if(OUTPUT == Output::modes::all || OUTPUT >= Output::modes::particles ){
                Output::particlesOutput(world, species, 1000);
            }
            if(world.getWallTime()>60){
                std::cout << "Time taken so far " << world.getWallTime()/60 << " minutes." << std::endl;
            }else{
                std::cout << "Time taken so far " << world.getWallTime() << " seconds." << std::endl;
            }
        }
        // std::cout << std::flush <<"\n";

       
    }
    if(world.getWallTime()>60){
        std::cout << "Simulation took " << world.getWallTime()/60 << " minutes." << std::endl;
    }else{
        std::cout << "Simulation took " << world.getWallTime() << " seconds." << std::endl;
    }
    return 0;
}
#endif

#ifdef SPHERE

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
        int solver_max_it = parseArgument(args, "--s_max_it", 2000); //or 2e4;
        type_calc phi_sphere = parseArgument(args, "--sphere_phi", -100.0);


        
        /*Instantiate World*/ 
        world_ptr = std::make_unique<World>(41, 21, 41, type_calc3{-0.2, -0.1, 0.0}, type_calc3{0.2, 0.1, 0.4});
        // world_ptr = std::make_unique<World>(41, 21, 41, type_calc3{-0.2,-0.1,0}, type_calc3{0.2,0.1,0.4});
        world_ptr ->setTimeStart();
        world_ptr->setTime(1e-7, 600);

        /*Instantiate sphere*/ 
        world_ptr->addObject<Sphere>(type_calc3(0, 0, 0.15), phi_sphere, 0.05);
        world_ptr->addObject<Rectangle>(type_calc3(0, 0, 0.35), phi_sphere, type_calc3(0.05, 0.05, 0.05));
        world_ptr->computeObjectID();
        std::string inlet_Face = "-z";
        world_ptr->addInlet(inlet_Face);

        // /*Instantiate species*/ 
        // species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 1e2);  // ions //last is mpw0
        // species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e1);
        // species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 2e3);

        // species.emplace_back("O+", 16*Const::amu, Const::q_e, *world_ptr, 1e12);  // ions //last is mpw0
        // species.emplace_back("O++", 16*Const::amu, 2*Const::q_e, *world_ptr, 5e11);
        species.emplace_back("O", 16*Const::amu, 0, *world_ptr, 1e12);
        
        // /*spherium (for rectangle also name temporary)*/
        // species.emplace_back("Sph", 100*Const::amu, 0, *world_ptr, 2e2);
        // species.emplace_back("Sph+", 100*Const::amu, Const::q_e, *world_ptr, 2e2);
        // species.emplace_back("e-", Const::m_e, -Const::q_e, *world_ptr);      // electrons


        /* Instantiate sources */
        const type_calc num_den_neutrals = 2e20; 
        const type_calc num_den_ions = 1e10; //mean ion density
        // const type_calc num_den_neutrals = 1e13; 
        // const type_calc num_den_ions = 1e10; //mean ion density
        type_calc T = 300;
        // sources.emplace_back(species[0], *world_ptr, 7000, 0.8*num_den_ions, T, inlet_Face); //O+
        // // sources.emplace_back(species[1], *world_ptr, 7000, 0.1*num_den_ions, T, inlet_Face); //O++
        // sources.emplace_back(std::make_unique<ColdBeamSource>(species[0], *world_ptr, 7000, 0.8*num_den_ions, inlet_Face)); //O+
        // sources.emplace_back(std::make_unique<ColdBeamSource>(species[1], *world_ptr, 7000, 0.1*num_den_ions, inlet_Face)); //O++
        // sources.emplace_back(std::make_unique<WarmBeamSource>(species[0], *world_ptr, 7000, 0.8*num_den_ions, T, inlet_Face)); //O+
        // sources.emplace_back(std::make_unique<WarmBeamSource>(species[1], *world_ptr, 7000, 0.1*num_den_ions, T, inlet_Face)); //O++
        sources.emplace_back(std::make_unique<WarmBeamSource>(species[0], *world_ptr, 7000, num_den_neutrals, 0, inlet_Face)); //O
        
        /* Instantiate interactions*/
        type_calc rate = 1e-4;
        // interactions.emplace_back(std::make_unique<ChemistryIonize>(species[3], species[4], *world_ptr, rate));
        interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));


        /*Instantiate solver*/ 
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        // solver_ptr->setReferenceValues(0, num_den_ions, 1.5);
        solver_ptr->setReferenceValues(0, 0, num_den_ions);

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
    Output::fieldsOutput(*world_ptr, species);

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
            sp.advanceSputtering(neutrals, spherium, world_ptr->getDt()); //TODO correct it so for species neutrals refer to species neutrals
            sp.computeNumberDensity();
            sp.sampleMoments();
            sp.computeMacroParticlesCount();
        }

        world_ptr->computeChargeDensity(species);
        // solver_ptr->solve();
        // solver_ptr->computeEF();

        if(world_ptr->getTs()>5){
            for(Species&  sp: species){
                sp.updateAverages();
            }
        }
        if(world_ptr->steadyState(species)){
            std::cout << "steady_state\n";
        }
        
        Output::diagOutput(*world_ptr, species);
        Output::screenOutput(*world_ptr, species);
        int ts = world_ptr->getTs();
        if(ts%10 == 0 || world_ptr->isLastTimeStep()){ //|| (ts > 140 && ts < 160)){
            Output::fieldsOutput(*world_ptr, species);
            Output::particlesOutput(*world_ptr, species, 10000);
            std::cout << "Time taken so far: " << world_ptr->getWallTime() << std::endl;

        }
       
    }
    std::cout << "Simulation took " << world_ptr->getWallTime() << " seconds." << std::endl;
    return 0;
}
#endif

#ifdef DISCHARGE
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
        int solver_max_it = parseArgument(args, "--s_max_it", 3800); //or 2e4;
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
        // ZROBIC dodać więcej czynników do samplowania prędkości i sprawdzić czy to pomaga
        species[2].loadParticleBoxThermal(type_calc3{world_ptr->getXc()[0], world_ptr->getXc()[1], world_ptr->getXc()[2] - 0.5*world_ptr->getL()[2]*0.85}, type_calc3(world_ptr->getL()[0]*0.01, world_ptr->getL()[1]*0.01, world_ptr->getL()[2]*0.02), num_den_electrons, T_eles);

        species[0].loadParticleBoxThermal(world_ptr->getXc(), type_calc3(world_ptr->getL()[0], world_ptr->getL()[1], world_ptr->getL()[2]*0.9), num_den_neutrals, T);
        
        /* Instantiate interactions*/
        interactions.emplace_back(std::make_unique<MC_MEX_Ionization>(species[0], species[1], species[2], *world_ptr));
        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[2], *world_ptr));
        //interactions.emplace_back(std::make_unique<DSMC_MEX>(species[1], *world_ptr));

        //interactions.emplace_back(std::make_unique<DSMC_MEX_Ionization>(species[1], species[2], species[3], *world_ptr));

        // type_calc rate_0_1 = 1e-2;
        // type_calc rate_1_2 = 1e-3;
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[0], species[1], species[3], *world_ptr, rate_0_1));
        // interactions.emplace_back(std::make_unique<ChemistryIonizeElectrons>(species[1], species[2], species[3], *world_ptr, rate_1_2));

        // interactions.emplace_back(std::make_unique<DSMC_MEX>(species[0], *world_ptr));


        /*Instantiate solver*/ 
        solver_ptr = std::make_unique<PotentialSolver>(*world_ptr, solver_max_it ,1e-4, static_cast<SolverType>(solver_type));  // Jeśli solver nie został zainicjalizowany
        //solver_ptr->setReferenceValues(0, num_den_ions, T); //simulationg electrons via boltzman relationship
        solver_ptr->setReferenceValues(0, 0, 1e20); //simulationg electrons directly 0 in ne zeroes infuence from boltzman relationship


        std::cout << "\nOptions:\n";
        std::cout << "Sybcycling: " << SUBCYCLING << "\n";
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
                    sp.sampleMoments();
                    sp.computeNumberDensity();
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
        //solver_ptr->solve();
        //solver_ptr->computeEF();


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

#endif
