#ifndef INSTANTIATOR_H
#define INSTANTIATOR_H
#include <memory>
#include "World.h"
#include "PotentialSolver.h"
#include "Species.h"
#include "funkc.h"



class Instantiator{
private:
    //world 
    int3                     Wnn    = {21, 21, 21};
    type_calc3               Wx0    = {-0.1, -0.1, 0.0};
    type_calc3               Wxm    = {0.1, 0.1, 0.2};
    type_calc                Wdt    = 2e-10;
    int                      Wts    = 1e4;
    std::vector<std::string> Wnames = {"grid", "x0", "xm", "dt", "ts"};

    //species and particles
    std::vector<std::string> Pids;
    std::vector<type_calc> Pmasses;
    std::vector<type_calc> Pcharges;
    std::vector<int3> Pgrids;
    std::vector<type_calc3> Pstarts;
    std::vector<type_calc3> Pends;
    std::vector<type_calc> Pdensities;
    bool QS = true;
    std::vector<std::string> Pnames = {"n", "m", "cha", "g", "x0", "xm", "d"};
    type_calc TEMP_MPW0 = 0;

    //solver
    std::vector<std::string> Snames     = {"mi", "tol"};
    SolverType Stype = PCG;

public:
    //solver 
    /*public bcoz Solver needs a reference to world constructed earlier, so we cant use anything else than constructor with
      reference to world object0, i might have confused reference with something else but point still stands */
    int       Smax_it    = 1e4;
    type_calc Stolerance = 1e-4;
    
    std::unique_ptr<World> world; //needen in readFunctions as a values storage xD 
    //std::vector<Species> species; /not needed
    //std::unique_ptr<PotentialSolver> solver; //not needed

    Instantiator(std::unique_ptr<World>& otherWorld, std::vector<Species>& otherSpecies, std::unique_ptr<PotentialSolver>& otherSolver);
    void readParametersWorld();
    void readParametersSpecies();
    void readParametersSolver();
};


#endif