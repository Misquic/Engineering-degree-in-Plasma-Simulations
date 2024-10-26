#ifndef POTENTIALSOLVER_H
#define POTENTIALSOLVER_H
//#include <>
#include "World.h"

class PotentialSolver{
protected:
    World&     world;         //reference to world object
    unsigned  max_solver_it;  //maximum number of solver iterations
    type_calc tolerance;      //tolerance used for convergence check

    type_calc phi0 = 0;
    type_calc n0 = 0;
    type_calc Te0 = 1;
public:
    /*constructors*/
    PotentialSolver(World& world, unsigned max_solver_it, type_calc tolerance);

    /*methods*/
    bool solveGS(); //solves potential using Gauss-Seidel and Succesive Over Relaxation: lokalna nadrelaksacja?
    void computeEF(); // computes electric field = -gradieng(phi)
    void setReferenceValues(type_calc phi0, type_calc n0, type_calc Te0);

};

#endif