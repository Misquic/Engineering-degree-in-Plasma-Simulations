#include "Interactions.h"

ChemistryIonize::ChemistryIonize(Species& neutrals, Species& ions, World& world, type_calc rate) noexcept: neutrals{neutrals},
 ions{ions}, world{world}, rate{rate}, dNi_field{world.ni, world.nj, world.nk}{
};

/*extended version with killing neutrals*/
void ChemistryIonize::apply(type_calc dt) noexcept{ //no electrons species //recombination is assumed to only occur on the surface when collision with object
    type_calc dV = world.getCellVolume();
    dNi_field.clear();
    int created = 0;
    /*loop over all cells*/
    for (int i=0; i<world.ni-1; i++){
		for (int j=0; j<world.nj-1; j++){ 
			for (int k=0; k<world.nk-1; k++){
                /*compute neutral and electron density at cell center*/
                type_calc den_a = neutrals.den.gather({i+0.5, j+0.5, k+0.5}); //number density of neutrals
                type_calc den_e = world.rho.gather({i+0.5, j+0.5, k+0.5})/Const::q_e; //number density of electrons
                type_calc d_den_i = rate*den_a*den_e*dt;

                dNi_field[i][j][k] = d_den_i; // <-MOJE

                int num_p = (int)(d_den_i*dV/ions.mpw0 + rnd()); //number of macroparticles to create
                if( num_p > den_a*dV/neutrals.mpw0){ //limitter
                    num_p = den_a*dV/neutrals.mpw0;
                }
                created+= num_p;
                for (int p = 0; p < num_p; p++){ //poprawic, zeby nie powstawalo wiecej niz jest
                    /*sample random particle in cell*/
                    type_calc3 lc{i+rnd(), j+rnd(), k+rnd()};
                    type_calc3 pos = world.LtoX(lc);
                    type_calc T = neutrals.T.gather(lc); //TODO check if its reliable
                    type_calc3 vel = neutrals.sampleV3th(T);
                    ions.addParticle(pos, vel, ions.mpw0);
                }
            }
        }
    }
    std::cout << "\nions created:" << created;
    ions.computeNumberDensity();


    // if(created > 0){
    //     std::cout << dNi_field << "\n";
    // }

    /*delete according number of particles*/ //poprawic, zeby dla kazdej celi usunac losowych tyle ile powstalo 
    std::vector<Particle>& particles = neutrals.getPartRef();
    int np = particles.size();
    int deleted = 0;
    type_calc P_sum = 0;
    for(int p = 0; p < np; p ++){
        Particle& part = particles[p];
        type_calc3 lc = world.XtoL(part.pos);
        type_calc dNa = dNi_field.gather(lc);

        int P = dNa*dV/neutrals.mpw0 + rnd();
        if(P>1){//kill
            deleted++;
            particles[p] = std::move(particles[np-1]);
            np--;
            p--;
            continue;
        }

    }
    std::cout << " neutrals deleted:" << deleted << "\n";

};

// /*booklike*/
// void ChemistryIonize::apply(type_calc dt) noexcept{ //no electrons species //recombination is assumed to only occur on the surface when collision with object
//     type_calc dV = world.getCellVolume();

//     /*loop over all cells*/
//     for (int i=0; i<world.ni-1; i++){
// 		for (int j=0; j<world.nj-1; j++){ 
// 			for (int k=0; k<world.nk-1; k++){
//                 /*compute neutral and electron density at cell center*/
//                 type_calc Na = neutrals.den.gather({i+0.5, j+0.5, k+0.5}); //number density of neutrals
//                 type_calc Ne = world.rho.gather({i+0.5, j+0.5, k+0.5})/Const::q_e; //number density of electrons
//                 type_calc dNi = rate*Na*Ne*dt;

//                 /*number of macroparticles to create*/
//                 int num_p = (int)(dNi*dV/ions.mpw0 + rnd());
//                 for (int p = 0; p < num_p; p++){
//                     /*sample random particle in cell*/
//                     type_calc3 lc{i+rnd(), j+rnd(), k+rnd()};
//                     type_calc3 pos = world.LtoX(lc);
//                     type_calc T = neutrals.T.gather(lc); //TODO check if its reliable
//                     type_calc3 vel = neutrals.sampleV3th(T);
//                     ions.addParticle(pos, vel, ions.mpw0);
//                 }
//             }
//         }
//     }
// };

////////////////////////////////// DSMC /////////////////////////////////////////

DSMC_MEX::DSMC_MEX(Species& species, World& world) noexcept: species{species}, world{world}{ //FOR INTERACTION BETWEEN ONE SPECIES

        m_reduced = species.mass*species.mass/(species.mass + species.mass); //TODO SEE BIRD'S APPENDIX A
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};

type_calc DSMC_MEX::evaluateSigma(type_calc v_rel) {
    return Const::pi * c[0] * c[0] * pow(c[2] / (v_rel*v_rel), c[1]-0.5) / c[3];

}

void DSMC_MEX::apply(type_calc dt) noexcept{ //FOR INTERACTION BETWEEN ONE SPECIES
    std::vector<std::vector<Particle*>> parts_in_cell =  species.sort_pointers(); // array of array of particle pointers

    type_calc dv = world.getCellVolume();
    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = species.mpw0;

    int n_collisions = 0;

    for(std::vector<Particle*>& parts: parts_in_cell){
        int np = parts.size();
        if(np <2) continue;

        type_calc n_groups_frac = 0.5 * np * np * mpw0 * sigma_v_rel_max * dt / dv;
        int n_groups = (int)(n_groups_frac + 0.5);

        for(int g = 0; g < n_groups; g++){
            int p1 = (int)(rnd()*np);
            int p2 = (int)(rnd()*np);
            while(p1 == p2){ //select random index of particle in cell
                p2 = (int)(rnd()*np);
            }

            type_calc v_rel = (parts[p1]->vel - parts[p2]->vel).length();

            type_calc sigma_v_rel = evaluateSigma(v_rel) * v_rel;
            if(sigma_v_rel > sigma_v_rel_max_temp) sigma_v_rel_max_temp = sigma_v_rel;

            type_calc P = sigma_v_rel/sigma_v_rel_max;

            if(P > rnd()){
                n_collisions++;
                collide(parts[p1]->vel, parts[p2]->vel, species.mass, species.mass);
            }
        }
    }
    std::cout << "ncollisions: " << n_collisions << "\n";
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
////////////////////////////////////////////////////////////////////////////////////////////

    // 	/*first we need to sort particles to cells*/
	// std::vector<Particle*> *parts_in_cell;
	// int n_cells = (world.ni-1)*(world.nj-1)*(world.nk-1);
	// parts_in_cell = new std::vector<Particle*> [n_cells];

	// /*sort particles to cells*/
	// for (Particle &part: species.getPartRef())
	// {
	// 	int c = world.XtoC(part.pos);
	// 	parts_in_cell[c].push_back(&part);
	// }

	// double sigma_cr_max_temp = 0;	/*reset for max computation*/
	// double dV = world.getCellVolume();	/*internal cell volume*/
	// double Fn = species.mpw0;	/*specific weight, using Bird's notation*/
	// int num_cols=0;	/*reset collision counter*/
					
	// /*now perform collisions*/
	// for (int c=0;c<n_cells;c++)
	// {
	// 	std::vector<Particle*> &parts = parts_in_cell[c];
	// 	int np = parts.size();
	// 	if (np<2) continue;

	// 	/*compute number of groups according to NTC*/
	// 	double ng_f = 0.5*np*np*Fn*sigma_v_rel_max*dt/dV;

	// 	int ng = (int)(ng_f+0.5);	/*number of groups, round*/
    //     //std::cout << "ng_f: " << ng_f << "\n";
	// 	/*assumes at least two particles per cell*/
	// 	for (int g=0;g<ng;g++)
	// 	{
	// 		int p1, p2;
	// 		p1 = (int)(rnd()*np);		/*returns some number between 0 and np-1 inclusive*/
		
	// 		do {
	// 			p2 = (int)(rnd()*np);
	// 		} while (p2==p1);

	// 		/*compute relative velocity*/
	// 		double3 cr_vec = parts[p1]->vel - parts[p2]->vel;
	// 		double cr = cr_vec.length();

	// 		/*evaluate cross section*/
	// 		double sigma = evaluateSigma(cr);

	// 		/*eval sigma_cr*/
	// 		double sigma_cr=sigma*cr;

	// 		/*update sigma_cr_max*/
	// 		if (sigma_cr>sigma_cr_max_temp)
	// 			sigma_cr_max_temp=sigma_cr;

	// 		/*eval prob*/
	// 		double P=sigma_cr/sigma_v_rel_max;

	// 		/*did the collision occur?*/
	// 		if (P>rnd())
	// 		{
	// 			num_cols++;
	// 			collide(parts[p1]->vel,parts[p2]->vel,species.mass, species.mass);
	// 		}
	// 	}
	// }

    // std::cout << "ncollisions: " << num_cols << "\n";
	// if (num_cols){
	// 	sigma_v_rel_max = sigma_cr_max_temp;
	// }



    // std::vector<std::vector<int>> parts_in_cell =  species.sort_indexes(); //array of array of particle indexes
}

void DSMC_MEX::collide(type_calc3& vel1, type_calc3& vel2, type_calc mass1, type_calc mass2){ // TODO implement better collisionn, not random isotropic
    type_calc3 cm = (mass1 * vel1 + mass2 * vel2)/(mass1 + mass2);

    type_calc3 v_rel = vel1 - vel2;
    type_calc v_rel_mag = v_rel.length();

    type_calc cos_ksi = 2*rnd() - 1;
    type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
    type_calc eps = 2*Const::pi * rnd();

    /*rotation*/
    v_rel[0] = v_rel_mag * cos_ksi;
    v_rel[1] = v_rel_mag * sin_ksi*std::cos(eps);
    v_rel[2] = v_rel_mag * sin_ksi*std::sin(eps);

    vel1 = cm + mass2/(mass1 + mass2)*v_rel;
    vel2 = cm - mass1/(mass1 + mass2)*v_rel;
};


