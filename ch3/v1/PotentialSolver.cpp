#include "PotentialSolver.h"

/*constructors*/
PotentialSolver::PotentialSolver(World& world, unsigned max_solver_it, type_calc tolerance): world(world), max_solver_it{max_solver_it}, tolerance{tolerance} {
};

/*methods*/
bool PotentialSolver::solveGS(){
    Field<type_calc>& phi = world.phi;  //references to avoid using world.phi
    Field<type_calc>& rho = world.rho;

    /*precalculate inverses*/
    static type_calc3 dx                = world.getDx();
    static type_calc  inv_d2x           = 1.0/(dx[0] * dx[0]);
    static type_calc  inv_d2y           = 1.0/(dx[1] * dx[1]);
    static type_calc  inv_d2z           = 1.0/(dx[2] * dx[2]);
    static type_calc  inv_eps_0         = 1.0/Const::eps_0;
    static type_calc  twos_over_invs    = 2.0*(inv_d2x + inv_d2y + inv_d2z);  // 2/d2x + 2/d2y + 2/d2z
    static type_calc  inv_twos_over_inv = 1.0/twos_over_invs;               // because multyplying is faster than deviding;

    /*for GS SOR*/
    type_calc new_phi = 0;                   
    const     type_calc SOR_weight = 1.4;  //weight used in SOR;

    /*for convergence checks*/
    type_calc L2{};                      //???
    type_calc sum{}; //
    type_calc R{}; //
    bool      converged = false;  

    unsigned it{};
    double ne{};
    /*solve potential using GS and SOR with neuman conditions*/
    /*GS*/


    for(it = 0; it < max_solver_it; it++){

        for(int i = 0; i < world.ni; i++){
            for(int j = 0; j < world.nj; j++){
                for(int k = 0; k < world.nk; k++){

                    if(world.object_id[i][j][k] > 0)
                        continue;

                    if(i == 0){
                        phi[i][j][k] = phi[i+1][j][k];
                    }                          // \ 
                          // |
                    else if(i == world.ni-1)            // |
                        phi[i][j][k] = phi[i-1][j][k];  // |
                    else if(j == 0)                     // |
                        phi[i][j][k] = phi[i][j+1][k];  // |
                    else if(j == world.nj-1)            // }- boundary
                        phi[i][j][k] = phi[i][j-1][k];  // |
                    else if(k == 0)                     // |
                        phi[i][j][k] = phi[i][j][k+1];  // |
                    else if(k == world.nk-1)            // |
                        phi[i][j][k] = phi[i][j][k-1];  // /
                    else{
                        //inside
                        //Boltzman relationship for electron density
                        ne = n0 * exp((phi[i][j][k] - phi0)/Te0);
                        new_phi = ((rho[i][j][k] - Const::q_e*ne)*inv_eps_0 + 
                        (phi[i-1][j][k]+phi[i+1][j][k])*inv_d2x + 
                        (phi[i][j-1][k]+phi[i][j+1][k])*inv_d2y + 
                        (phi[i][j][k-1]+phi[i][j][k+1])*inv_d2z) * inv_twos_over_inv;
                        //if(it < 2 && k == 1) std::cerr << "phi[i-1]: " << phi[i-1][j][k] << " phi[i+1]: "<< phi[i+1][j][k] << " " << (phi[i-1][j][k]+phi[i+1][j][k])*inv_d2x << " ";
                        /*SOR*/
                        phi[i][j][k] = phi[i][j][k] + SOR_weight * (new_phi - phi[i][j][k]);
                    }

                }   
            }
        }
        /*check for convergence*/
        if(it%25 == 0){
            //std::cerr << phi[0][1][9] << std::endl;
            sum = 0;
            for(int i = 1; i < world.ni-1; i++){
                for(int j = 1; j < world.nj-1; j++){
                    for(int k = 1; k < world.nk-1; k++){
                        if(world.object_id[i][j][k] > 0) continue;
                        if(i == 0)
                            R = phi[i][j][k] - phi[i+1][j][k];
                        else if(i == world.ni-1)
                            R = phi[i][j][k] - phi[i-1][j][k];
                        else if(j == 0)
                            R = phi[i][j][k] - phi[i][j+1][k];
                        else if(j == world.nj-1)
                            R = phi[i][j][k] - phi[i][j-1][k];
                        else if(k == 0)
                            R = phi[i][j][k] - phi[i][j][k+1];
                        else if(k == world.nk-1)
                            R = phi[i][j][k] - phi[i][j][k-1];
                        else{
                            R = -phi[i][j][k]*twos_over_invs + 
                            (rho[i][j][k] - Const::q_e*ne)*inv_eps_0 + 
                            (phi[i-1][j][k]+phi[i+1][j][k])*inv_d2x + 
                            (phi[i][j-1][k]+phi[i][j+1][k])*inv_d2y + 
                            (phi[i][j][k-1]+phi[i][j][k+1])*inv_d2z;
                        }
                        sum += R*R;
                    }
                }
            }
            L2 = sqrt(sum / (world.ni * world.nj * world.nk));
            if (L2 < tolerance){
                converged = true;
                break;
            }
        }
    }

    if(!converged){
        std::cerr << "GS SOR failed to converge, L2 = " << L2 << " tolerance = " << tolerance << " time step = " << world.getTs() << std::endl;
    }

    //std::cerr << "Iterations made: " << it << " L2 = " << L2 << std::endl; 

    return converged;
};
void PotentialSolver::computeEF(){
    Field<type_calc>& phi = world.phi;  //references to avoid using world.phi
    Field<type_calc3>& ef = world.ef;  //references to avoid using world.phi
    type_calc3* ef_ptr = nullptr; 

    /*precalculate inverses*/
    type_calc3 dx                = world.getDx();
    type_calc  inv_2dx           = 1.0/(2 * dx[0]);
    type_calc  inv_2dy           = 1.0/(2 * dx[1]);
    type_calc  inv_2dz           = 1.0/(2 * dx[2]);

    for(int i = 0; i < world.ni; i++){ //?write explicite without if? for performance
        for(int j = 0; j < world.nj; j++){
            for(int k = 0; k < world.nk; k++){
                ef_ptr = &ef[i][j][k]; //pointer to vec3 at i j k used kinda like reference in solve, but to avoid allocating in each loop

                if(i == 0){ //Ex
                    //finite difference forward, second order
                    (*ef_ptr)[0] = (3*phi[i][j][k] - 4*phi[i+1][j][k] + phi[i+2][j][k]) * inv_2dx;
                }
                else if(i == world.ni-1 )
                {
                    //finite difference backward, second order
                    (*ef_ptr)[0] = (-phi[i-2][j][k] + 4*phi[i-1][j][k] - 3*phi[i][j][k]) * inv_2dx;
                }
                else{
                    //finite difference central
                    (*ef_ptr)[0] = (phi[i-1][j][k] - phi[i+1][j][k]) * inv_2dx;
                };
                
                if(j == 0){ //Ey
                    //finite difference forward, second order
                    (*ef_ptr)[1] = (3*phi[i][j][k] - 4*phi[i][j+1][k] + phi[i][j+2][k]) * inv_2dy;
                }
                else if(j == world.nj-1 )
                {
                    //finite difference backward, second order
                    (*ef_ptr)[1] = (-phi[i][j-2][k] + 4*phi[i][j-1][k] - 3*phi[i][j][k]) * inv_2dy;
                }
                else{
                    //finite difference central
                    (*ef_ptr)[1] = (phi[i][j-1][k] - phi[i][j+1][k]) * inv_2dy;
                };

                if(k == 0){ //Ez
                    //finite difference forward, second order
                    (*ef_ptr)[2] = (3*phi[i][j][k] - 4*phi[i][j][k+1] + phi[i][j][k+2]) * inv_2dz;
                }
                else if(k == world.nk-1 )
                {
                    //finite difference backward, second order
                    (*ef_ptr)[2] = (-phi[i][j][k-2] + 4*phi[i][j][k-1] - 3*phi[i][j][k]) * inv_2dz;
                }
                else{
                    //finite difference central
                    (*ef_ptr)[2] = (phi[i][j][k-1] - phi[i][j][k+1]) * inv_2dz;
                };
            }
        }
    }
}

void PotentialSolver::setReferenceValues(type_calc phi0, type_calc n0, type_calc Te0){
    this->phi0 = phi0;
    this->n0 = n0;
    this->Te0 = Te0;
};
