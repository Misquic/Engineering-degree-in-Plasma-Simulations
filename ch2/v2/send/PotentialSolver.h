#ifndef POTENTIALSOLVER_H
#define POTENTIALSOLVER_H
//#include <>
#include "World.h"

class PotentialSolver{
protected:
    World&     world;         //reference to world object
    unsigned  max_solver_it;  //maximum number of solver iterations
    type_calc tolerance;      //tolerance used for convergence check

public:
    /*constructors*/
    PotentialSolver(World& world, unsigned max_solver_it, type_calc tolerance);

    /*methods*/
    bool solveGS(); //solves potential using Gauss-Seidel and Succesive Over Relaxation: lokalna nadrelaksacja?
    void computeEF(); // computes electric field = -gradieng(phi)


};

#endif