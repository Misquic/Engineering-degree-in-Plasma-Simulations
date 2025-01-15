#ifndef POTENTIALSOLVER_H
#define POTENTIALSOLVER_H
//#include <>
#include "World.h"
#include "Matrix.h"

namespace vec {
    tcvector deflate(const Field<type_calc>& f3); //convert 3D Field to 1D vector
    void inflate(const tcvector& v, Field<type_calc>& f3); //convert 1D vector to 3D Field
    template<class data_type>
    void inflate(const std::vector<data_type>& v, Field<data_type>& f3){
    for( int i = 0; i < f3.ni; i++){
        for( int j = 0; j < f3.nj; j++){
            for(int k = 0; k < f3.nk; k++){
                f3[i][j][k] = v[k*f3.ni*f3.nj + j*f3.ni + i];
            }
        }
    }
};
    type_calc dot(const tcvector& v, const tcvector& u); //dot product in higer dimension
    type_calc norm(const tcvector& v); //norm or length in higher dimension
}

enum SolverType{GS, PCG, QN};

class PotentialSolver{
protected:
    World&     world;         //reference to world object
    const SolverType solver_type;   //type of used solver

    /*GS*/
    unsigned  PCG_max_solver_it;  //maximum number of solver iterations
    unsigned  GS_max_solver_it;  //maximum number of solver iterations
    type_calc tolerance;      //tolerance used for convergence check
    type_calc phi0 = 0;       //reference phi for boltzman relationship
    type_calc n0 = 0;         //reference number of micropart for boltzman relationship
    type_calc Te0 = 1;        //reference electron temperature for boltzman relationship
    const type_calc SOR_weight = 1.4;
    
    /*NRPCG*/
    Matrix A;

    enum             NodeType{REGULAR, NEUMANN, DIRICHLET};  //enum for node types
    std::vector<int> node_type;                              //flag for different node types
public:
    /*constructors*/
    PotentialSolver(World& world, unsigned max_solver_it, type_calc tolerance, SolverType solver_type);

    /*methods*/
    bool solve(); 
    bool solveGS();                                                        //solves potential using Gauss-Seidel and Succesive Over Relaxation
    
    bool solveNRPCG();                                                     //solves potential using Newton-Rapshon Preconditioned Conjugate Gradient uses PCGlinear and GSlinear if it fails to converge
    bool solvePCGlinear(const Matrix &A, tcvector &x, const tcvector &b);  // solves matrix equation Ax = b using Preconditioned Conjugate Gradient
    bool solveGSlinear(const Matrix& A, tcvector& x, const tcvector& b);   //solves matrix equation Ax = b Gauss-Seidel and Succesive Over Relaxation
   
    bool solveQN();                                                        //solves quasi neutral system

    void computeEF();                                                      // computes electric field = -gradieng(phi)
    void setReferenceValues(type_calc phi0, type_calc n0, type_calc Te0);

    unsigned get_GS_max_it();
    unsigned get_PCG_max_it();


private:
    void buildMatrix();
    void precalculate();
    /*precalculate values*/
    /*potential*/
    type_calc3 dx;
    type_calc  inv_d2x;
    type_calc  inv_d2y;
    type_calc  inv_d2z;
    type_calc  inv_eps_0;
    type_calc  twos_over_inv_d2;
    type_calc  inv_twos_over_inv;
    /*electic field*/
    type_calc  inv_2dx;
    type_calc  inv_2dy;
    type_calc  inv_2dz;
    /*build matrix*/
    type_calc  inv_dx;
    type_calc  inv_dy;
    type_calc  inv_dz;
    /*PCG*/
    Matrix M; //inverse of Jacobi matrix
};

std::ostream& operator<<(std::ostream& out, SolverType& type);
std::istream& operator>>(std::istream& in, SolverType& type);

#endif