#include "PotentialSolver.h"
#include <cmath>
#include "funkc.h"
#ifdef DEBUG
#include "Outputs.h"
#endif
/////////////////////////////////////// namespace vec /////////////////////////////////////////////////////

tcvector vec::deflate(const Field<type_calc>& f3){
    tcvector v(f3.ni * f3.nj * f3.nk);
    for( int i = 0; i < f3.ni; i++){
        for( int j = 0; j < f3.nj; j++){
            for(int k = 0; k < f3.nk; k++){
                v[f3.U(i,j,k)] = f3[i][j][k];
            }
        }
    }
    return v;
};
void vec::inflate(const tcvector& v, Field<type_calc>& f3){
    for( int i = 0; i < f3.ni; i++){
        for( int j = 0; j < f3.nj; j++){
            for(int k = 0; k < f3.nk; k++){
                f3[i][j][k] = v[f3.U(i,j,k)];
            }
        }
    }
};
type_calc vec::dot(const tcvector& v, const tcvector& u){
    double dot = 0;
    size_t n_unknowns = v.size();
    for (size_t i = 0; i < n_unknowns; i++)
        dot += v[i] * u[i];
    return dot;
};
type_calc vec::norm(const tcvector& v){
    return sqrt(dot(v,v)/v.size());
}

/////////////////////////////////////// Potential Solver //////////////////////////////////////////////////


/*constructors*/
PotentialSolver::PotentialSolver(World& world, unsigned max_solver_it, type_calc tolerance, SolverType solver_type): world(world),
 PCG_max_solver_it{max_solver_it}, tolerance{tolerance}, A{world.ni*world.nj*world.nk}, M{world.ni*world.nj*world.nk}, solver_type{solver_type} {
    GS_max_solver_it = 20*max_solver_it;
    if(solver_type == GS){
        GS_max_solver_it = max_solver_it;
    }
    precalculate();
    buildMatrix();
};

/*methods*/
bool PotentialSolver::solve(){
    switch(solver_type)    {
    case GS:
        return solveGS();
    case PCG:
        return solveNRPCG();
    case QN:
        return solveQN();
    default:
        std::cerr << "Wrong solver type\n";
        return false;
    }
};

bool PotentialSolver::solveGS(){
    Field<type_calc>& phi = world.phi;  //references to avoid using world.phi
    Field<type_calc>& rho = world.rho;

    /*for GS SOR*/
    type_calc new_phi{};                   
    // const     type_calc SOR_weight = 1.4;  //weight used in SOR;

    /*for convergence checks*/
    type_calc L2{};                      
    type_calc sum{}; 
    type_calc R{}; 
    bool      converged = false;  

    unsigned it{};
    double ne{};
    /*solve potential using GS and SOR with neuman conditions*/
    /*GS*/
    for(it = 0; it < GS_max_solver_it; it++){

        for(int i = 0; i < world.ni; i++){
            for(int j = 0; j < world.nj; j++){
                for(int k = 0; k < world.nk; k++){

                    if(world.object_id[i][j][k] > 0)
                        continue;

                    if(i == 0)
                        phi[i][j][k] = phi[i+1][j][k];  // 
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
                        /*SOR*/
                        phi[i][j][k] = phi[i][j][k] + SOR_weight * (new_phi - phi[i][j][k]);
                    }

                }   
            }
        }
        /*check for convergence*/
        if(it%25 == 0){
            sum = 0;
            for(int i = 0; i < world.ni; i++){
                for(int j = 0; j < world.nj; j++){
                    for(int k = 0; k < world.nk; k++){
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
                            ne = n0*exp((phi[i][j][k] - phi0)/Te0);
                            R = -phi[i][j][k]*twos_over_inv_d2 + 
                            (rho[i][j][k] - Const::q_e*ne)*inv_eps_0 + 
                            (phi[i-1][j][k]+phi[i+1][j][k])*inv_d2x + 
                            (phi[i][j-1][k]+phi[i][j+1][k])*inv_d2y + 
                            (phi[i][j][k-1]+phi[i][j][k+1])*inv_d2z;
                        }
                        sum += R*R;
                    }
                }
            }
            L2 = std::sqrt(sum / (world.nv));
            if (L2 < tolerance){
                converged = true;
                break;
            }
        } // it%25
    }

    if(!converged){
        std::cerr << "GS SOR failed to converge, L2 = " << L2 << " tolerance = " << tolerance << " time step = " << world.getTs() << std::endl;
    }
    return converged;
};
//??? pointers slower slightly due to dereference?
bool checkVec(const tcvector& vec){
    bool result = false;
    for(int i = 0; i < vec.size(); i++){
        if(std::isnan(vec[i])){
            std::cerr << "nan in vec:" << i << "\n";
            result = true;
        }
    }
    return result;
}
bool PotentialSolver::solveNRPCG(){
    const int NR_MAX_IT = 20; //maximum number of NR iterations
    const type_calc NR_TOL = 1e-3;
    int n_unknowns = A.n_unknowns;

    Matrix   J(n_unknowns);  //Jacobian
    tcvector P(n_unknowns);  //Jacobi Preconditioner
    tcvector y(n_unknowns);  //solution Jy = F
    tcvector phi_ = vec::deflate(world.phi);
    tcvector rho_ = vec::deflate(world.rho);

    for(int i = 0; i < n_unknowns; i++){
        if(node_type[i] == NEUMANN) rho_[i] = 0;
        else if(node_type[i] == DIRICHLET) rho_[i] = phi_[i];
        else rho_[i] = -rho_[i]*inv_eps_0;
    }

    type_calc norm{};
    bool converged = false;

    for(int it = 0; it < NR_MAX_IT; it++){
        std::cout << "\rNRPCG starting iteration nr" << it+1 << "\n";
        tcvector F = A*phi_-rho_; //calculate F by first subtracting linear term

        for(int i = 0; i < n_unknowns; i++){ // subtract b(x) on regular nodes
            if(node_type[i] == REGULAR){
                F[i] -= Const::q_e*n0*exp((phi_[i]-phi0)/Te0)*inv_eps_0;
            }
        }

        for(int i = 0; i < n_unknowns; i++){
            if(node_type[i] == REGULAR){
                P[i] = Const::q_e*n0/(Const::eps_0*Te0)*exp((phi_[i]-phi0)/Te0); //TODO merge with above?
                // if(std::isnan(Const::q_e*n0/(Const::eps_0*Te0)*exp((phi_[i]-phi0)/Te0))){
                //     std::cerr << "nan";
                // }
            }
        }

        J = A.diagSubtract(P);
        #ifdef DEBUG
            Output::convergenceOutput(it, world.getTs());
        
        #endif
        if(!solvePCGlinear(J,y,F)){
            solveGSlinear(J,y,F);
        }
        for(int i = 0; i < n_unknowns; i++){
            if(node_type[i] == DIRICHLET) y[i] = 0;
        }
        phi_ = phi_-y; // TODO use =-
        norm = vec::norm(y);
        if(norm<NR_TOL){
            converged = true;
            break;
        }
    }

    if(!converged){
        std::cerr << "NR+PCG failed to converge, norm = " << norm << " tolerance = " << NR_TOL << " time step = " << world.getTs() << std::endl;
    }
    vec::inflate(phi_, world.phi);
    return converged;
};

bool PotentialSolver::solvePCGlinear(const Matrix& A, tcvector& x, const tcvector& b){ // sometimes doesn't converge, if that happens it makes solution worse, sometimes GS cannot make it right fast enough then it can lead to nans
    bool converged = false;
    type_calc L2{};
    //Matrix M = A.invDiagonal();

    // checkVec(x);
    // checkVec(b);
    tcvector g = A*x-b;
    // checkVec(g);
    tcvector s = M*g;
    // checkVec(s);
    tcvector d = -1*s;
    // checkVec(d);

    //preallocate
    tcvector z;z.reserve(d.size());

    // const unsigned logInterval = std::max(1u, PCG_max_solver_it / 100);

    for(unsigned int it = 0; it < PCG_max_solver_it; it ++){
#ifdef DEBUG
        if(it%int(PCG_max_solver_it*0.01)==0){
        }
#endif
        z = A*d;                                                                          //matrix mult
        // checkVec(d);
        // checkVec(z);
        type_calc alpha = g*s;
        type_calc beta = d*z;
        if(beta == 0){
            std::cerr << "PotentialSolver::solvePCGlinear: beta = 0" << beta << "\n";
            throw std::runtime_error("division by 0, beta = 0\n");
            
        }
#ifdef DEBUG
#endif

        x = x + (alpha/beta)*d;
        // checkVec(x);
        g = g + (alpha/beta)*z;
        // checkVec(g);
        s = M*g;                                                                          //matrix mult
        // checkVec(s);
        beta = alpha;
        alpha = g*s;
        d = (alpha/beta)*d - s;
        // checkVec(d);
        L2 = vec::norm(g);
        // if(it%logInterval==0){
        //     std::cout << "\r                                      \r"<< type_calc(it)/PCG_max_solver_it*100 << "% PCG, L2: " << L2;
        // }    

#ifdef DEBUG
        Output::convergenceOutput(L2,it, world.getTs());
#endif

        if(L2 < tolerance){
            converged = true;
            break;
        }
    }    

    std::cout << "\r                                      \r"<< 100 << "% PCG" << L2 << "\n";
    if(!converged){
        std::cerr << "PCGlinear failed to converge, norm(g) = " << L2 << " tolerance = " << tolerance << " time step = " << world.getTs()<< std::endl;

    }
    return converged;

};

bool PotentialSolver::solveGSlinear(const Matrix &A, tcvector& x, const tcvector& b){ //always makes solution better
    type_calc L2{};
    bool converged = false; //sprawdzic x czy jest nan 

    for(int it = 0; it < GS_max_solver_it; it++){
        for(int i = 0; i < A.n_unknowns; i++){

            type_calc S = A.multiplyRow(i,x) - A(i,i)*x[i];
            type_calc phi_new = (b[i] - S)/A(i,i);
            
            x[i] = x[i] + 1.0 * (phi_new - x[i]) ;//change 1.0 to SOR_WEIGHT? 
           
        }
        if(it%25 == 0){
            tcvector R = A*x-b;
            L2 = vec::norm(R);
            #ifdef DEBUG
            if(it%int(GS_max_solver_it*0.01)==0){
                std::cout << "\r                                      \r"<< type_calc(it)/GS_max_solver_it*100 << "% GS, L2: " << L2;;
            }
            #endif

            if(L2<tolerance){
                converged = true;
                break;
            }
        }
    }
    if(!converged){
        std::cerr <<  "GSlinear failed to converge, L2 = " << L2 << " tolerance = " << tolerance << " time step = " << world.getTs()<< std::endl;
    }
    return converged;
}

bool PotentialSolver::solveQN(){
    std::cerr << "quasi neutral not implemented yet\n";
    return false;
};

void PotentialSolver::computeEF(){
    Field<type_calc>& phi = world.phi;
    Field<type_calc3>& ef = world.ef;
    type_calc3* ef_ptr = nullptr; 

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
};
void PotentialSolver::setReferenceValues(type_calc phi0, type_calc n0, type_calc Te0){
    this->phi0 = phi0;
    this->n0 = n0;
    this->Te0 = Te0;
};
void PotentialSolver::buildMatrix(){
    node_type = std::vector<int>(A.n_unknowns);
    

    for(int k = 0; k < world.nk; k++){ //why indexex running in reverse? //cash locality perhaps
        for(int j = 0; j < world.nj; j++){
            for(int i = 0; i < world.ni; i++){
                int u = k*world.nj*world.ni + j*world.ni + i;
                A.clearRow(u);

                if(world.object_id[i][j][k] > 0 || world.node_type[i][j][k] == DIRICHLET){
                    A(u,u) = 1;
                    node_type[u] = DIRICHLET;
                    continue;
                }
                
                node_type[u] = NEUMANN;
                if(i == 0){
                    A(u,u) = inv_dx;
                    A(u,u+1) = -inv_dx;
                }
                else if(i == world.ni-1){
                    A(u,u) = inv_dx;
                    A(u,u-1) = -inv_dx;
                }
                else if(j == 0){
                    A(u,u) = inv_dy;
                    A(u,u+world.ni) = -inv_dy;
                }
                else if(j == world.nj -1){
                    A(u,u) = inv_dy;
                    A(u,u-world.ni) = -inv_dy;
                }
                else if(k == 0){
                    A(u,u) = inv_dz;
                    A(u,u+world.nj*world.ni) = -inv_dz;
                }
                else if(k == world.nk -1){
                    A(u,u) = inv_dz;
                    A(u,u-world.nj*world.ni) = -inv_dz;
                }
                else{
                    node_type[u] = REGULAR;
                    A(u,u-world.ni*world.nj) = inv_d2x;
                    A(u,u-world.ni) = inv_d2y;
                    A(u,u-1) = inv_d2z;
                    A(u,u) = -twos_over_inv_d2;
                    A(u,u+1) = inv_d2z;
                    A(u,u+world.ni) = inv_d2y;
                    A(u,u+world.ni*world.nj) = inv_d2x;
                }
            }
        }
    }
    // precalculate M, inverse of Jacobi matrix
    vec::inflate(node_type, world.node_type);
    M = A.invDiagonal();
};

void PotentialSolver::precalculate(){
    /*potential*/
    dx = world.getDx();
    inv_d2x = 1.0/(dx[0] * dx[0]);
    inv_d2y = 1.0/(dx[1] * dx[1]);
    inv_d2z = 1.0/(dx[2] * dx[2]);
    inv_eps_0 = 1.0/Const::eps_0;
    twos_over_inv_d2 = 2.0*(inv_d2x + inv_d2y + inv_d2z);
    inv_twos_over_inv = 1.0/twos_over_inv_d2;
    /*electric field*/
    inv_2dx = 1.0/(2 * dx[0]);
    inv_2dy = 1.0/(2 * dx[1]);
    inv_2dz = 1.0/(2 * dx[2]);
    /*build matrix*/
    inv_dx = 1.0/dx[0];
    inv_dy = 1.0/dx[1];
    inv_dz = 1.0/dx[2];
    /*PCG*/
}

std::ostream& operator<<(std::ostream& out, SolverType& type){
    switch(type){
    case SolverType::GS:
        out << "GS";
        break;
    case SolverType::PCG:
        out << "PCG";
        break;
    case SolverType::QN:
        out << "QN";
        break;
    default:
        break;
    }
    return out;
};

std::istream& operator>>(std::istream& in, SolverType& type){
    std::string input;
    in >> input;

    if (input == "GS") {
        type = SolverType::GS;
    } else if (input == "PCG") {
        type = SolverType::PCG;
    } else if (input == "QN") {
        type = SolverType::QN;
    } else {
        throw(std::invalid_argument("Wrong Solver Type, setting SolverType::GS"));
        type = GS;  // DEFAULT
    }

    return in;

};

unsigned PotentialSolver::get_GS_max_it(){
    return GS_max_solver_it;
};
unsigned PotentialSolver::get_PCG_max_it(){
    return PCG_max_solver_it;
};

