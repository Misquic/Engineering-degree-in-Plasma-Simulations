#include "PotentialSolver.h"
#include <cmath>
#include "funkc.h"

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
    size_t nu = v.size();
    for (size_t j=0;j<nu;j++)
        dot+=v[j]*u[j];
    return dot;
    //return v*u;
};
type_calc vec::norm(const tcvector& v){
    double sum = 0;
    int nu = v.size();
    for (int j=0;j<nu;j++)
        sum+=v[j]*v[j];
    return sqrt(sum/nu);
    //return sqrt((v*v)/v.size());
}


/////////////////////////////////////// Potential Solver //////////////////////////////////////////////////


/*constructors*/
PotentialSolver::PotentialSolver(World& world, unsigned max_solver_it, type_calc tolerance, SolverType solver_type): world(world),
 max_solver_it{max_solver_it}, tolerance{tolerance}, A{world.ni*world.nj*world.nk}, M{world.ni*world.nj*world.nk}, solver_type{solver_type} {
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

    /*precalculate inverses*/ //TODO add as world or solver parameters -> speed!?
    // static type_calc3 dx                = world.getDx();
    // static type_calc  inv_d2x           = 1.0/(dx[0] * dx[0]);
    // static type_calc  inv_d2y           = 1.0/(dx[1] * dx[1]);
    // static type_calc  inv_d2z           = 1.0/(dx[2] * dx[2]);
    // static type_calc  inv_eps_0         = 1.0/Const::eps_0;
    // static type_calc  twos_over_invs    = 2.0*(inv_d2x + inv_d2y + inv_d2z);  // 2/d2x + 2/d2y + 2/d2z
    // static type_calc  inv_twos_over_inv = 1.0/twos_over_invs;               // because multyplying is faster than deviding;

    /*for GS SOR*/
    type_calc new_phi{};                   
    // const     type_calc SOR_weight = 1.4;  //weight used in SOR;

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

    //std::cerr << "Iterations made: " << it << " L2 = " << L2 << std::endl; 

    return converged;
};

bool PotentialSolver::solveNRPCG(){
    // const int NR_MAX_IT = 20;
    // const type_calc NR_TOL = 1e-3;
    // int n_unknowns = A.n_unknowns;

    // Matrix J(n_unknowns);
    // tcvector Q(n_unknowns);
    // tcvector y(n_unknowns);
    // tcvector x = vec::deflate(world.phi);
    // tcvector b = vec::deflate(world.rho);

    // for(int i = 0; i < n_unknowns; i++){
    //     if(node_type[i] == NEUMANN) b[i] = 0;
    //     else if(node_type[i] == DIRICHLET) b[i] = x[i];
    //     else b[i] = -b[i]*inv_eps_0;
    // }

    // type_calc norm{};
    // bool converged = false;

    // for(int it = 0; it < NR_MAX_IT; it++){
    //     tcvector F = A*x-b;

    //     for(int i = 0; i < n_unknowns; i++){
    //         if(node_type[i] == REGULAR){
    //             F[i] -= Const::q_e*n0*exp((x[i]-phi0)/Te0)*inv_eps_0;
    //         }
    //     }
    //     for(int i = 0; i < n_unknowns; i++){
    //         if(node_type[i] == REGULAR){
    //             Q[i] = Const::q_e*n0/(Const::eps_0*Te0)*exp((x[i]-phi0)/Te0); //TODO merge with above?
    //         }
    //     }

    //     J = A.diagSubtract(Q);
    //     if(!solvePCGlinear(J,y,F)){
    //         solveGSlinear(J,y,F);
    //     }
    //     for(int i = 0; i < n_unknowns; i++){
    //         if(node_type[i] == DIRICHLET) y[i] = 0;
    //     }
    //     x = x-y;
    //     norm = vec::norm(y);
    //     if(norm<NR_TOL){
    //         converged = true;
    //         break;
    //     }
    // }

    // if(!converged){
    //     std::cerr << "NR+PCG failed to converge, norm = " << norm << " tolerance = " << NR_TOL << " time step = " << world.getTs() << std::endl;
    // }
    // //std::cerr << "Iterations made: " << it << " norm = " << norm << std::endl; 
    // vec::inflate(x, world.phi);
    // return converged;
    /*main NR iteration loop*/
	const int NR_MAX_IT=20;		/*maximum number of NR iterations*/
	const double NR_TOL = 1e-3;
	int nu = A.n_unknowns;

	Matrix J(nu);
	tcvector P(nu);
	tcvector y(nu);
	tcvector x = vec::deflate(world.phi);
	tcvector b = vec::deflate(world.rho);

	/*set RHS to zero on boundary nodes (zero electric field)
      and to existing potential on fixed nodes */
    for (int u=0;u<nu;u++)
    {
		if (node_type[u]==NEUMANN) b[u] = 0;			/*neumann boundary*/
        else if (node_type[u]==DIRICHLET) b[u] = x[u];	/*dirichlet boundary*/
        else b[u] = -b[u]/Const::eps_0;            /*regular node*/
    }

	double norm;
	bool converged=false;
	for(int it=0;it<NR_MAX_IT;it++)
	{
		/*compute F by first subtracting the linear term */
		tcvector F = A*x-b;

		/*subtract b(x) on regular nodes*/
		for (int n=0;n<nu;n++)
			if (node_type[n]==REGULAR)	/*regular nodes*/
				F[n] -= Const::q_e*n0*exp((x[n]-phi0)/Te0)/Const::eps_0;

		/*Compute P, diagonal of d(bx)/dphi*/
		for (int n=0;n<nu;n++)
		{
			if (node_type[n]==REGULAR)
				P[n] = n0*Const::q_e/(Const::eps_0*Te0)*exp((x[n]-phi0)/Te0);
		}

		/*Compute J = A-diag(P)*/
		Matrix J = A.diagSubtract(P);

		/*solve Jy=F*/
		if (!solvePCGlinear(J,y,F))
			solveGSlinear(J,y,F);

		/*clear any numerical noise on Dirichlet nodes*/
		for (int u=0;u<nu;u++)
			if (node_type[u]==DIRICHLET) y[u]=0;

		/*x=x-y*/
		x = x-y;

		norm=vec::norm(y);
		//cout<<"NR norm: "<<norm<<endl;

		if (norm<NR_TOL)
		{
			converged=true;
			break;
		}
	}

	if (!converged)
		std::cout<<"NR+PCG failed to converge, norm = "<<norm<<std::endl;

	/*convert to 3d data*/
	vec::inflate(x,world.phi);
	return converged;
};
bool PotentialSolver::solvePCGlinear(const Matrix& A, tcvector& x, const tcvector& b){
    // bool converged = false;
    // tcvector tmp = A*x; 
    // tcvector g = A*x-b; //
    // for(int i = 0; i < A.n_unknowns; i ++){
    //     if(std::isnan(tmp[i]) || std::isnan(x[i]) || std::isnan(b[i])){
    //         std::cout << " " << (A*x)[i] << " " << x[i] << " " << b[i] << " " << g[i];
    //         std::cout  << " " << i << " " << world.getTs() << "\n";
    //         throw(std::invalid_argument("PCG problem"));
    //     }
    // }
    // tcvector s = M*g; //
    // tcvector d = -1*s; //
    // type_calc L2{};
    // //Matrix M = A.invDiagonal();

    // for(unsigned it = 0; it < max_solver_it; it ++){
    //     tcvector z = A*d;//
    //     type_calc alpha = g*s;
    //     type_calc beta = d*z; //


    //     x = x + (alpha/beta)*d;
    //     g = g + (alpha/beta)*z;
    //     s = M*g;

    //     beta = alpha;
    //     alpha = g*s;
    //     // tu gdzieś dochodzi do nan
    //     d = (alpha/beta)*d - s;
    //     L2 = vec::norm(g);
        
    //     if(L2 < tolerance){
    //         converged = true;
    //         break;
    //     }
    // }    

    // if(!converged){
    //     std::cerr << "PCGlinear failed to converge, norm(g) = " << L2 << " tolerance = " << tolerance << " time step = " << world.getTs() << std::endl;
    //     //std::cout  << "\n" << " g: " << g << std::endl;
    // }
    // //std::cerr << "Iterations made: " << it << " L2 = " << L2 << std::endl; 
    // return converged;

    bool converged= false;

	double l2 = 0;
	Matrix M = A.invDiagonal(); //inverse of Jacobi preconditioner

	/*initialization*/
	tcvector g = A*x-b;
	tcvector s = M*g;
	tcvector d = -1*s;

	for (unsigned it=0;it<max_solver_it;it++)
	{
		tcvector z = A*d;
		double alpha = vec::dot(g,s);
		double beta = vec::dot(d,z);

		x = x+(alpha/beta)*d;
		g = g+(alpha/beta)*z;
		s = M*g;

		beta = alpha;
		alpha = vec::dot(g,s);

		d = (alpha/beta)*d-s;
		l2 = vec::norm(g);
		if (l2<tolerance) {converged=true;break;}
	}

	if (!converged)	std::cerr<<"PCG failed to converge, norm(g) = "<<l2<<std::endl;
    return converged;

}
//////////////////// AAAAAAAAAAAAAAAAAAAAAAAA //////////////////////
//////////////////// AAAAAAAAAAAAAAAAAAAAAAAA //////////////////////
bool PotentialSolver::solveGSlinear(const Matrix &A, tcvector& x, const tcvector& b){
    type_calc L2{};
    bool converged = false; //sprawdzic x czy jest nan 

    for(int it = 0; it < max_solver_it; it++){
        for(int i = 0; i < A.n_unknowns; i++){

            type_calc S = A.multiplyRow(i,x) - A(i,i)*x[i]; //tu coś nie tak?
            type_calc phi_new = (b[i] - S)/A(i,i);
            
            x[i] = x[i] + 1.0 * (phi_new - x[i]) ;//change 1.0 to SOR_WEIGHT? //czy tu sie robi nan?
           
        }
        if(it%25 == 0){
            tcvector R = A*x-b;
            L2 = vec::norm(R);
            if(L2<tolerance){
                converged = true;
                break;
            }
        }
    }
    if(!converged){
        std::cerr <<  "GSlinear failed to converge, L2 = " << L2 << " tolerance = " << tolerance << " time step = " << world.getTs() << std::endl;
    }
    return converged;
}

bool PotentialSolver::solveQN(){
    std::cerr << "quasi neutral not implemented yet\n";
    return false;
};

void PotentialSolver::computeEF(){
    Field<type_calc>& phi = world.phi;  //references to avoid using world.phi
    Field<type_calc3>& ef = world.ef;  //references to avoid using world.phi
    type_calc3* ef_ptr = nullptr; 

    // /*precalculate inverses*/
    // type_calc3 dx                = world.getDx();
    // type_calc  inv_2dx           = 1.0/(2 * dx[0]);
    // type_calc  inv_2dy           = 1.0/(2 * dx[1]);
    // type_calc  inv_2dz           = 1.0/(2 * dx[2]);

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
    // node_type = std::vector<int>(A.n_unknowns);
    // //node_type.reserve(A.n_unknowns);

    // for(int k = 0; k < world.nk; k++){ //why indexex running in reverse?
    //     for(int j = 0; j < world.nj; j++){
    //         for(int i = 0; i < world.ni; i++){
    //             int u = k*world.nj*world.ni + j*world.ni + i;
    //             A.clearRow(u);

    //             if(world.object_id[i][j][k] > 0){
    //                 A(u,u) = 1;
    //                 node_type[u] = DIRICHLET;
    //                 continue;
    //             }
                
    //             node_type[u] = NEUMANN;
    //             if(i == 0){
    //                 A(u,u) = inv_dx;
    //                 A(u,u+1) = -inv_dx;
    //             }
    //             else if(i == world.ni-1){
    //                 A(u,u) = inv_dx;
    //                 A(u,u-1) = -inv_dx;
    //             }
    //             else if(j == 0){
    //                 A(u,u) = inv_dy;
    //                 A(u,u+world.ni) = -inv_dy;
    //             }
    //             else if(j == world.nj -1){
    //                 A(u,u) = inv_dy;
    //                 A(u,u-world.ni) = -inv_dy;
    //             }
    //             else if(k == 0){
    //                 A(u,u) = inv_dz;
    //                 A(u,u+world.nj*world.ni) = -inv_dz;
    //             }
    //             else if(k == world.nk -1){
    //                 A(u,u) = inv_dz;
    //                 A(u,u-world.nj*world.ni) = -inv_dz;
    //             }
    //             else{
    //                 node_type[u] = REGULAR;
    //                 A(u,u-world.ni*world.nj) = inv_d2x;
    //                 A(u,u-world.ni) = inv_d2y;
    //                 A(u,u-1) = inv_d2z;
    //                 A(u,u) = -twos_over_invs;
    //                 A(u,u+1) = inv_d2z;
    //                 A(u,u+world.ni) = inv_d2y;
    //                 A(u,u+world.ni*world.nj) = inv_d2x;
    //             }
    //         }
    //     }
    // }
    // // precalculate M, inverse of Jacobi matrix
    // vec::inflate(node_type, world.node_type);
    // M = A.invDiagonal();
    double3 dh = world.getDx();
	double idx = 1.0/dh[0];
	double idy = 1.0/dh[1];
	double idz = 1.0/dh[2];
    double idx2 = idx*idx;	/*1/(dx*dx)*/
	double idy2 = idy*idy;
	double idz2 = idz*idz;
	int ni = world.ni;
	int nj = world.nj;
	int nk = world.nk;
	int nu = ni*nj*nk;

	/*reserve space for node types*/
	node_type.reserve(nu);

	/*solve potential*/
	for (int k=0;k<nk;k++)
        for (int j=0;j<nj;j++)
        	for (int i=0;i<ni;i++)
            {
                int u = world.phi.U(i,j,k);
                A.clearRow(u);
                //dirichlet node?
				if (world.object_id[i][j][k]>0)
                {
                    A(u,u)=1;	//set 1 on the diagonal
                    node_type[u] = DIRICHLET;
                    continue;
                }

				//Neumann boundaries
				node_type[u] = NEUMANN;		//set default
                if (i==0) {A(u,u)=idx;A(u,u+1)=-idx;}
                else if (i==ni-1) {A(u,u)=idx;A(u,u-1)=-idx;}
                else if (j==0) {A(u,u)=idy;A(u,u+ni)=-idy;}
                else if (j==nj-1) {A(u,u)=idy;A(u,u-ni)=-idy;}
                else if (k==0) {A(u,u)=idz;A(u,u+ni*nj)=-idz;}
				else if (k==nk-1) {
					A(u,u)=idz;
					A(u,u-ni*nj)=-idz;}
                else {
                	//standard internal stencil
                	A(u,u-ni*nj) = idz2;
                	A(u,u-ni) = idy2;
                	A(u,u-1) = idx2;
                	A(u,u) = -2.0*(idx2+idy2+idz2);
                	A(u,u+1) = idx2;
                	A(u,u+ni) = idy2;
                	A(u,u+ni*nj) = idz2;
                	node_type[u] = REGULAR;	//regular internal node
                }
            }
	//solveQN();
};

void PotentialSolver::precalculate(){
    /*potential*/
    dx = world.getDx();
    inv_d2x = 1.0/(dx[0] * dx[0]);
    inv_d2y = 1.0/(dx[1] * dx[1]);
    inv_d2z = 1.0/(dx[2] * dx[2]);
    inv_eps_0 = 1.0/Const::eps_0;
    twos_over_invs = 2.0*(inv_d2x + inv_d2y + inv_d2z);
    inv_twos_over_inv = 1.0/twos_over_invs;
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

