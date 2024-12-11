#include <algorithm>
#include "Interactions.h"

ChemistryIonize::ChemistryIonize(Species& neutrals, Species& ions, World& world, type_calc rate) noexcept: neutrals{neutrals},
 ions{ions}, world{world}, rate{rate}, dNi_field{world.ni, world.nj, world.nk} {
};
void ChemistryIonize::apply(type_calc dt) noexcept{ //no electrons species //recombination is assumed to only occur on the surface when collision with object
    type_calc dV = world.getCellVolume();
    int created = 0;
    std::vector<std::vector<Particle*>> sorted_neutrals =  neutrals.sort_pointers(); // array of array of particle pointers
    std::vector<Particle*> to_remove;
    std::vector<Particle> ions_particles = ions.getPartRef();

    int ni = world.ni;
    int nj = world.nj;
    /*loop over all cells*/
    for (int i=0; i<world.ni-1; i++){
		for (int j=0; j<world.nj-1; j++){ 
			for (int k=0; k<world.nk-1; k++){
                /*compute neutral and electron density at cell center*/
                type_calc den_a = neutrals.den.gather({i+0.5, j+0.5, k+0.5}); //number density of neutrals
                type_calc den_e = world.rho.gather({i+0.5, j+0.5, k+0.5})/Const::q_e; //number density of electrons
                type_calc d_den_i = rate*den_a*den_e*dt;

                int c_index = k*(ni-1)*(nj-1) + j*(ni-1) + i; //cell index
                int num_part_in_cell = sorted_neutrals[c_index].size();

                //std::cout << "num to create: " << d_den_i*dV/ions.mpw0 << "\n";
                int num_p_to_create = (int)(d_den_i*dV/ions.mpw0 + rnd()); //number of macroparticles to create
                if( num_p_to_create > num_part_in_cell){ //limitter
                    num_p_to_create = num_part_in_cell;
                }
                created+= num_p_to_create;
                to_remove.reserve(to_remove.size() + num_p_to_create);
                int np = neutrals.getConstPartRef().size(); // ??
                if(num_p_to_create > 0){
                    ions_particles.reserve(ions_particles.size() + num_p_to_create);
                }

                for (int p = 0; p < num_p_to_create; p++){
                    
                    /*sample random particle in cell*/
                    int pic_index = rnd(0,num_part_in_cell);
                    Particle* samped_neutral_ptr = sorted_neutrals[c_index][pic_index];
                    
                    /*create ion and electron*/
                    type_calc3 lc = world.XtoL(samped_neutral_ptr->pos);
                    type_calc T = neutrals.T.gather(lc);
                    ions.addParticle(samped_neutral_ptr->pos, samped_neutral_ptr->vel, ions.mpw0);
                    
                    /*delete sampled neutral*/
                    to_remove.push_back(samped_neutral_ptr);
                    sorted_neutrals[c_index][pic_index] = std::move(sorted_neutrals[c_index][num_part_in_cell-1]);
                    num_part_in_cell--; // CHECK IF IT doesnt step over 0
                }
            }
        }
    }
    auto& particles = neutrals.getPartRef();
    particles.erase( std::remove_if( particles.begin(), particles.end(), [&to_remove](const Particle& part) { return std::any_of( to_remove.begin(), to_remove.end(), [&part](Particle* ptr) { return &part == ptr; }); }),particles.end());
    
    std::cout << "\nions created:" << created;
};


ChemistryIonizeElectrons::ChemistryIonizeElectrons(Species& neutrals, Species& ions, Species& electrons, World& world, type_calc rate) noexcept: neutrals{neutrals},
 ions{ions}, electrons{electrons}, world{world}, rate{rate}, dNi_field{world.ni, world.nj, world.nk} {
};


/*extended version with killing neutrals*/
void ChemistryIonizeElectrons::apply(type_calc dt) noexcept{ //no electrons species //recombination is assumed to only occur on the surface when collision with object
    type_calc dV = world.getCellVolume();
    //dNi_field.clear();
    int created = 0;
    std::vector<std::vector<Particle*>> sorted_neutrals =  neutrals.sort_pointers(); // array of array of particle pointers
    std::vector<Particle*> to_remove;
    std::vector<Particle>& ions_particles = ions.getPartRef();
    std::vector<Particle>& electrons_particles = electrons.getPartRef();

    int ni = world.ni;
    int nj = world.nj;
    /*loop over all cells*/
    for (int i=0; i<world.ni-1; i++){
		for (int j=0; j<world.nj-1; j++){ 
			for (int k=0; k<world.nk-1; k++){
                /*compute neutral and electron density at cell center*/
                type_calc den_a = neutrals.den.gather({i+0.5, j+0.5, k+0.5}); //number density of neutrals
                // type_calc den_e = world.rho.gather({i+0.5, j+0.5, k+0.5})/Const::q_e; //number density of electrons
                type_calc den_e = electrons.den.gather({i+0.5, j+0.5, k+0.5}); //number density of electrons
                type_calc d_den_i = rate*den_a*den_e*dt;

                int c_index = k*(ni-1)*(nj-1) + j*(ni-1) + i; //cell index
                int num_part_in_cell = sorted_neutrals[c_index].size();

                //std::cout << "num to create: " << d_den_i*dV/ions.mpw0 << "\n";
                int num_p_to_create = (int)(d_den_i*dV/ions.mpw0 + rnd()); //number of macroparticles to create
                if( num_p_to_create > num_part_in_cell){ //limitter
                    num_p_to_create = num_part_in_cell;
                }
                created+= num_p_to_create;
                to_remove.reserve(to_remove.size() + num_p_to_create);
                int np = neutrals.getConstPartRef().size(); // ??
                if(num_p_to_create > 0){
                    ions_particles.reserve(ions_particles.size() + num_p_to_create);
                    electrons_particles.reserve(electrons_particles.size() + num_p_to_create);
                }


                for (int p = 0; p < num_p_to_create; p++){
                    
                    /*sample random particle in cell*/
                    int pic_index = rnd(0,num_part_in_cell);
                    Particle* samped_neutral_ptr = sorted_neutrals[c_index][pic_index];
                    
                    /*create ion and electron*/
                    type_calc3 lc = world.XtoL(samped_neutral_ptr->pos);
                    type_calc T = neutrals.T.gather(lc);
                    ions.addParticle(samped_neutral_ptr->pos, samped_neutral_ptr->vel, ions.mpw0);
                    electrons.addParticle(samped_neutral_ptr->pos, samped_neutral_ptr->vel, electrons.mpw0); // change it to random v3th?
                    
                    /*delete sampled neutral*/
                    to_remove.push_back(samped_neutral_ptr); // CHANGE it to mpw marking
                    sorted_neutrals[c_index][pic_index] = std::move(sorted_neutrals[c_index][num_part_in_cell-1]);
                    num_part_in_cell--; // CHECK IF IT doesnt step over 0
                }
            }
        }
    }
    auto& particles = neutrals.getPartRef();
    particles.erase( std::remove_if( particles.begin(), particles.end(), [&to_remove](const Particle& part) { return std::any_of( to_remove.begin(), to_remove.end(), [&part](Particle* ptr) { return &part == ptr; }); }),particles.end());
    
    std::cout << "\nions created:" << created;
};

////////////////////////////////// DSMC /////////////////////////////////////////

DSMC_MEX::DSMC_MEX(Species& species, World& world) noexcept: species1{species}, species2{species}, world{world}{ //FOR INTERACTION BETWEEN ONE SPECIES
        apply_strategy = [this](type_calc dt){applyOneSpecies(dt);};


        m_reduced = species1.mass*species2.mass/(species1.mass + species2.mass); //TODO SEE BIRD'S APPENDIX A
        dv = world.getCellVolume();
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};
DSMC_MEX::DSMC_MEX(Species& species1, Species& species2, World& world): species1{species1}, species2{species2}, world{world}{ //FOR INTERACTION BETWEEN TWO SPECIES
        if(species1.mpw0 != species2.mpw0){
            throw("species must have the same macroparticle weight for this algorithm to work properly.");
        }
        apply_strategy = [this](type_calc dt){applyTwoSpecies(dt);};


        m_reduced = species1.mass*species2.mass/(species1.mass + species2.mass); //TODO SEE BIRD'S APPENDIX A
        dv = world.getCellVolume();
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};

type_calc DSMC_MEX::evaluateSigma(type_calc v_rel) {
    return Const::pi * c[0] * c[0] * pow(c[2] / (v_rel*v_rel), c[1]-0.5) / c[3];

}

void DSMC_MEX::apply(type_calc dt) noexcept{
    apply_strategy(dt);
}

void DSMC_MEX::applyOneSpecies(type_calc dt) noexcept{ //FOR INTERACTION BETWEEN ONE SPECIES
    std::vector<std::vector<Particle*>> parts_in_cell =  species1.sort_pointers(); // array of array of particle pointers

    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = species1.mpw0;

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
                collide(parts[p1]->vel, parts[p2]->vel, species1.mass, species1.mass);
            }
        }
    }
    std::cout << "ncollisions: " << n_collisions << "\n";
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
}

void DSMC_MEX::applyTwoSpecies(type_calc dt) noexcept{ //FOR INTERACTION BETWEEN TWO SPECIES
    std::vector<std::vector<Particle*>> parts_in_cell1 =  species1.sort_pointers(); // array of array of particle pointers
    std::vector<std::vector<Particle*>> parts_in_cell2 =  species2.sort_pointers(); // array of array of particle pointers

    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = species1.mpw0;

    int n_collisions = 0;
    
    int num_cells = (world.ni-1) * (world.nj -1) * (world.nk-1);
    for(int c = 0; c < num_cells; c++){
        std::vector<Particle*>& parts1 = parts_in_cell1[c];
        std::vector<Particle*>& parts2 = parts_in_cell2[c];
        int np1 = parts1.size();
        int np2 = parts2.size();
        if(np1 < 1 || np2 < 1) continue;

        type_calc n_groups_frac = 0.5 * np1 * np2 * mpw0 * sigma_v_rel_max * dt / dv;
        int n_groups = (int)(n_groups_frac + 0.5);

        for(int g = 0; g < n_groups; g++){
            int p1 = (int)(rnd()*np1);
            int p2 = (int)(rnd()*np2);

            type_calc v_rel = (parts1[p1]->vel - parts2[p2]->vel).length();

            type_calc sigma_v_rel = evaluateSigma(v_rel) * v_rel;
            if(sigma_v_rel > sigma_v_rel_max_temp) sigma_v_rel_max_temp = sigma_v_rel;

            type_calc P = sigma_v_rel/sigma_v_rel_max;

            if(P > rnd()){
                n_collisions++;
                collide(parts1[p1]->vel, parts2[p2]->vel, species1.mass, species1.mass);
            }
        }
    }
    std::cout << "ncollisions: " << n_collisions << "\n";
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
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

////////////////////////////////////// DSMC_IONISATION ///////////////////////////////////////////
DSMC_MEX_IONIZATION::DSMC_MEX_IONIZATION(Species& neutrals, Species& ions, Species& electrons, World& world): neutrals{neutrals}, ions{ions},
 electrons{electrons}, world{world}{ //FOR INTERACTION BETWEEN ONE SPECIES
        if(neutrals.mpw0 != electrons.mpw0){
            throw std::invalid_argument("species must have the same macroparticle weight for this algorithm to work properly.");
        }
        if(!neutrals.E_ion || neutrals.E_ion < 0){
            throw std::invalid_argument("neutral species must have proper ionization energy E_ion (in Joules)");
        }
        E_ion = neutrals.E_ion; //It should be scaled by mpw0 right?
        dv = world.getCellVolume();
        m_reduced = neutrals.mass*electrons.mass/(neutrals.mass + electrons.mass); //TODO SEE BIRD'S APPENDIX A
        c[0] = 1.016e-18;
        c[1] = 9.824e+00;
        c[2] = 6.175e+01;
        c[3] = 666;

		// c[0] = 4.07e-10;
		// c[1] = 0.77;
		// c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		// c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};
void DSMC_MEX_IONIZATION::apply(type_calc dt) noexcept{
    std::vector<std::vector<std::unique_ptr<Particle>>> neutrals_in_cell = neutrals.sort_pointers_unique(); // array of array of particle pointers
    std::vector<std::vector<std::unique_ptr<Particle>>> electrons_in_cell = electrons.sort_pointers_unique(); // array of array of particle pointers


    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = neutrals.mpw0;

    int n_collisions = 0;
    int n_ionizations = 0;
    
    int num_cells = (world.ni-1) * (world.nj -1) * (world.nk-1);

    for(int c = 0; c < num_cells; c++){
        // std::vector<Particle*>& parts_neutrals = neutrals_in_cell[c];
        // std::vector<Particle*>& parts_electrons = electrons_in_cell[c];

        int np_neus = neutrals_in_cell[c].size();
        int np_eles = electrons_in_cell[c].size();
        if(np_neus < 1 || np_eles < 1) continue;

        type_calc n_groups_frac = 0.5 * np_neus * np_eles * mpw0 * sigma_v_rel_max * dt / dv;
        int n_groups = (int)(n_groups_frac + 0.5);
        if(n_groups > np_neus){
            n_groups = np_neus-1;
        }
        //std::cout << "n_groups_frac: " << n_groups_frac << " ";
        for(int g = 0; g < n_groups; g++){
            //std::cout << "np_neus: " << np_neus << " np_eles: " << np_eles << "\n";
            int p_neus = (int)(rnd()*np_neus);
            while(neutrals_in_cell[c][p_neus]->macro_weight==0){
                p_neus = (int)(rnd()*np_neus); 
            }

            int p_eles = (int)(rnd()*np_eles);

            type_calc v_rel = (neutrals_in_cell[c][p_neus]->vel - electrons_in_cell[c][p_eles]->vel).length();

            type_calc sigma_v_rel = evaluateSigma(v_rel) * v_rel;
            if(sigma_v_rel > sigma_v_rel_max_temp) sigma_v_rel_max_temp = sigma_v_rel;

            type_calc P = sigma_v_rel/sigma_v_rel_max;

            //std::cout << "P: " << P << " sigma_rel: " << sigma_v_rel << " n_groups_frac: " << n_groups_frac << " np_neus: " << np_neus << " np_eles: " << np_eles << " sigma_v_rel_max: "<< sigma_v_rel_max << "\n";

            bool ionised = false;
            if(P > rnd()){
                n_collisions++;
                ionised = collide(neutrals_in_cell[c][p_neus]->vel, electrons_in_cell[c][p_eles]->vel, neutrals.mass, electrons.mass);
                if(ionised){
                    n_ionizations++;

                    int ions_to_create = (int)(ions.mpw0/neutrals.mpw0 + rnd());
                    for(int i = 0; i < ions_to_create; i++){
                        ions.addParticle(neutrals_in_cell[c][p_neus]->pos, neutrals_in_cell[c][p_neus]->vel, ions.mpw0);
                    }

                    // type_calc3 v_rel = parts_neutrals[p_neus]->vel - parts_electrons[p_eles]->vel;
                    // type_calc v_rel_mag = v_rel.length();
                    // type_calc E_rel = 0.5 * neutrals.mass * v_rel_mag * v_rel_mag; //in Joules
                    //type_calc3 lc = world.XtoL(parts_neutrals[p_neus]->pos);
                    //type_calc T = electrons.T.gather(lc);
                    electrons.addParticle(neutrals_in_cell[c][p_neus]->pos, electrons.sampleV3th(100), electrons.mpw0); // change it to random v3th?

                    neutrals_in_cell[c][p_neus]->macro_weight = 0;
                    // if(p_neus != np_neus-1){
                    //     neutrals_in_cell[c][p_neus] = std::move(neutrals_in_cell[c][np_neus-1]);
                    // }
                    // np_neus--;
                }
            }
        }
    }
    std::vector<Particle>& particles = neutrals.getPartRef();
    int np = particles.size();
    for(int p = 0; p < np; p++){
        if(!particles[p].macro_weight){
            particles[p] = std::move(particles[np-1]);
            p--;
            np--;
        }
    }
    particles.erase(particles.begin() + np, particles.end());
    
    std::cerr << "n_collisions: " << n_collisions << " n_ionizations: " << n_ionizations << "\n";
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
};
type_calc DSMC_MEX_IONIZATION::evaluateSigma(type_calc v_rel){
    type_calc E = 0.5 * neutrals.mass * v_rel * v_rel / Const::q_e;
    // type_calc E = 0.5 * electrons.mass * v_rel * v_rel;
    return 1e14*c[0] * std::log(E/c[1])/E*std::exp(-c[2]/E); //fit to data obtained via LXcat Source: RBEQ Fit to [Brook et al. J. Phys. B: At. Mol. Phys. 11, 3115 (1978)], [Rothe et al. Phys. Rev. 125, 582 (1962)], [Zipf, Planet Space Sci 33, 1303 (1985)], [Thompson et al. J. Phys. B: At. Mol. Opt. Phys. 28, 1321 (1995)]. Updated: 3 January 2024.
};
bool DSMC_MEX_IONIZATION::collide(type_calc3& vel1, type_calc3& vel2, type_calc mass1, type_calc mass2){
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

    type_calc E_rel = 0.5 * mass1 * v_rel_mag * v_rel_mag; //in Joules
    //std::cout << " E_rel: " << E_rel << " E_ion: " << E_ion << "\n";
    // type_calc E_rel = 0.5 * mass2 * v_rel_mag * v_rel_mag; //in Joules // TRY!!
    type_calc P = E_rel/E_ion;
    if(E_rel > E_ion) return true;
    return false;
};

