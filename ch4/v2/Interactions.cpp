#include <algorithm>
#include <unordered_map>
#include <ranges>
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
        mass1 = species.mass;
        mass2 = species.mass;
        sum_mass = mass1 + mass2;
        dv = world.getCellVolume();
		c[0] = 4.07e-10;
		c[1] = 0.77;
		c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};
DSMC_MEX::DSMC_MEX(Species& species1, Species& species2, World& world): species1{species1}, species2{species2}, world{world}{ //FOR INTERACTION BETWEEN TWO SPECIES
        if(species1.mpw0 != species2.mpw0){
            throw std::invalid_argument("species must have the same macroparticle weight for this algorithm to work properly.");
        }
        apply_strategy = [this](type_calc dt){applyTwoSpecies(dt);};


        m_reduced = species1.mass*species2.mass/(species1.mass + species2.mass); //TODO SEE BIRD'S APPENDIX A
        mass1 = species1.mass;
        mass2 = species2.mass;
        sum_mass = mass1 + mass2;
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
                collide(parts[p1]->vel, parts[p2]->vel);
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
                collide(parts1[p1]->vel, parts2[p2]->vel);
            }
        }
    }
    std::cout << "ncollisions: " << n_collisions << "\n";
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
}

void DSMC_MEX::collide(type_calc3& vel1, type_calc3& vel2){ // TODO implement better collisionn, not random isotropic
    type_calc3 cm = (mass1 * vel1 + mass2 * vel2)/sum_mass;

    type_calc3 v_rel = vel1 - vel2;
    type_calc v_rel_mag = v_rel.length();

    type_calc cos_ksi = 2*rnd() - 1;
    type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
    type_calc eps = 2*Const::pi * rnd();

    /*rotation*/
    v_rel[0] = v_rel_mag * cos_ksi;
    v_rel[1] = v_rel_mag * sin_ksi*std::cos(eps);
    v_rel[2] = v_rel_mag * sin_ksi*std::sin(eps);

    vel1 = cm + mass2/sum_mass*v_rel;
    vel2 = cm - mass1/sum_mass*v_rel;
};

////////////////////////////////////// DSMC_IONISATION ///////////////////////////////////////////
// ZROBIC sprawdzić jednostki
// ZROBIC afscra folder szybszy
// czy jak zmienić gęstość neutralnych i elektronów  i rate poszukać z literatury (zacząć od blogu) w podejściu chemicznym to czy pomoże
// 
DSMC_MEX_IONIZATION::DSMC_MEX_IONIZATION(Species& neutrals, Species& ions, Species& electrons, World& world): neutrals{neutrals}, ions{ions},
 electrons{electrons}, world{world}{ //FOR INTERACTION BETWEEN NEUTRALS AND ELECTRONS
        if(neutrals.mpw0 != electrons.mpw0){
            throw std::invalid_argument("species must have the same macroparticle weight for this algorithm to work properly.");
        }
        if(!neutrals.E_ion || neutrals.E_ion < 0){
            throw std::invalid_argument("neutral species must have proper ionization energy E_ion (in Joules)");
        }
        if(neutrals.mpw0 < ions.mpw0){
            std::cerr << "neutrals.mpw0: " << neutrals.mpw0 << " ions.mpw0: " << ions.mpw0;
            throw std::invalid_argument("neutral species must have greater mpw0 than ions to prperly calculate number of ions to create during ionisation, best if neutrals.mpw0 = k*ions.mpw0, k in natural numbers");
        }
        E_ion_times_10 = 10*neutrals.E_ion; //It should be not scaled by mpw0 right?
        dv = world.getCellVolume();
        m_reduced = neutrals.mass*electrons.mass/(neutrals.mass + electrons.mass); //TODO SEE BIRD'S APPENDIX A
        mass_neus = neutrals.mass;
        mass_eles = electrons.mass;
        sum_mass = mass_neus + mass_eles;
        num_cells = (world.ni-1) * (world.nj -1) * (world.nk-1);

        c[0] = 1.016e-18;
        c[1] = 9.824e+00;
        c[2] = 6.175e+01;
        c[3] = 666; //not used

		// c[0] = 4.07e-10;
		// c[1] = 0.77;
		// c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		// c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};


void DSMC_MEX_IONIZATION::apply(type_calc dt) noexcept{ // test for no pointers
    // std::vector<std::vector<int>> neutrals_in_cell = neutrals.sort_indexes(); // array of array of particle pointers
    // std::vector<std::vector<int>> electrons_in_cell = electrons.sort_indexes(); // array of array of particle pointers
    std::unordered_map<int, std::vector<int>> neutrals_in_cell_index; //will be rare because indexes only for when there are both electrons and neutrons
    std::unordered_map<int, std::vector<int>> electrons_in_cell_index; //rare
    electrons.map_indexes(electrons_in_cell_index); // implement sorting neutrals in map but only for cells that already have electrons
    neutrals.map_indexes(neutrals_in_cell_index, electrons_in_cell_index); //sorts and makes it rarer because only when there are both cpecies in cell
    std::vector<Particle>& parts_neutrals = neutrals.getPartRef();
    std::vector<Particle>& parts_electrons = electrons.getPartRef();

    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = neutrals.mpw0;
    int n_collisions = 0;
    int n_ionizations = 0;
    #ifdef DEBUG
    std::cerr << "\n";
    #endif
    for(const std::pair<int, std::vector<int>>& pair: electrons_in_cell_index ){
        int c = pair.first;
        const std::vector<int>& indexes_electrons = pair.second;//electrons_in_cell_index[c];
    
        std::vector<int>& indexes_neutrals = neutrals_in_cell_index[c];
        int np_neus_in_cell = indexes_neutrals.size();
        int np_eles_in_cell = indexes_electrons.size();

        type_calc n_groups_frac = 0.5 * np_neus_in_cell * np_eles_in_cell * mpw0 * sigma_v_rel_max * dt / dv;
        #ifdef DEBUG
        std::cerr << "n_groups_frac: " << n_groups_frac << " c: " << c;
        #endif
        int n_groups = (int)(n_groups_frac + 0.5);
        if(n_groups > np_neus_in_cell){ //cannot ionize more then it is aviable to ionize
            n_groups = np_neus_in_cell-1;
        }
        //std::cout << "n_groups_frac: " << n_groups_frac << " ";
        for(int g = 0; g < n_groups; g++){
            //std::cout << "np_neus_in_cell: " << np_neus_in_cell << " np_eles_in_cell: " << np_eles_in_cell << "\n";
            // int p_neus = (int)(rnd()*np_neus_in_cell); //sample random neutral particle
            int indexes_neutral_sample_index = rnd(0,np_neus_in_cell);
            int p_neu = indexes_neutrals[indexes_neutral_sample_index]; //sample random neutral particle
            int checked_neutrals = 0;
            while(parts_neutrals[p_neu].macro_weight==0 && checked_neutrals < 300){ //if it is marked to kill
                p_neu = indexes_neutrals[rnd(0,np_neus_in_cell)]; //should always find some but maximum is 300
                checked_neutrals++;
            }
            if(checked_neutrals >= 300){
                std::cerr << "particle to collide not found in cell (exceeded 300 checks threshold)\n";
                break;
            }
            int p_ele = indexes_electrons[rnd(0,np_eles_in_cell)]; //sample random electron particle

            type_calc v_rel = (parts_neutrals[p_neu].vel - parts_electrons[p_ele].vel).length();
            type_calc sigma_v_rel = evaluateSigma(v_rel) * v_rel;
            if(sigma_v_rel > sigma_v_rel_max_temp)
                sigma_v_rel_max_temp = sigma_v_rel;
            type_calc P = sigma_v_rel/sigma_v_rel_max;
            #ifdef DEBUG
            std::cerr << " P: " << P << " sigma_v_rel: " << sigma_v_rel << " sigma_ve_rel_max: " << sigma_v_rel_max << "\n";
            #endif
            bool ionised = false;
                if(P > rnd()){
                n_collisions++;
                type_calc E_new_ele = 0;
                collide(parts_neutrals[p_neu].vel, parts_electrons[p_ele].vel, ionised, E_new_ele);
                if(ionised){
                    n_ionizations++;
                    int ions_to_create = (int)(neutrals.mpw0/ions.mpw0 + rnd()); //assuming neutrals.mpw0 >= ions.mpw0
                    for(int i = 0; i < ions_to_create; i++){
                        ions.addParticle(parts_neutrals[p_neu].pos, parts_neutrals[p_neu].vel, ions.mpw0);
                    }

                    // type_calc3 v_rel = parts_neutrals[p_neus]->vel - parts_electrons[p_eles]->vel;
                    // type_calc v_rel_mag = v_rel.length();
                    // type_calc E_rel = 0.5 * neutrals.mass * v_rel_mag * v_rel_mag; //in Joules
                    // type_calc3 lc = world.XtoL(parts_neutrals[p_neus]->pos);
                    // type_calc T = electrons.T.gather(lc);
                    type_calc T = E_new_ele/Const::k;
                    #ifdef DEBUG
                    std::cerr << "T new electron: " << T << "\n";
                    #endif
                    electrons.addParticle(parts_neutrals[p_neu].pos, electrons.sampleV3th(T), electrons.mpw0); // change it to random v3th?

                    parts_neutrals[p_neu].macro_weight = 0; //mark nuetral to kill

                    if(indexes_neutral_sample_index != np_neus_in_cell-1){
                        indexes_neutrals[indexes_neutral_sample_index] = std::move(indexes_neutrals[np_neus_in_cell-1]); //kill ionised neutral index in cell
                    }else{
                        indexes_neutrals.pop_back();
                    }
                    np_neus_in_cell--;
                }
            }
        }
    }
    if(n_collisions){//change it to acc es only ones to kill
        int np = parts_neutrals.size();
        for(int p = 0; p < np; p++){
            if(parts_neutrals[p].macro_weight == 0){
                parts_neutrals[p] = std::move(parts_neutrals[np-1]);
                p--;
                np--;
            }
        }
        parts_neutrals.erase(parts_neutrals.begin() + np, parts_neutrals.end());
    }
    
    std::cerr << "n_collisions: " << n_collisions << " n_ionizations: " << n_ionizations;
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
};

type_calc DSMC_MEX_IONIZATION::evaluateSigma(type_calc v_rel){
    type_calc E = 0.5 * m_reduced * v_rel * v_rel / Const::q_e;
    // type_calc E = 0.5 * electrons.mass * v_rel * v_rel;
    // 1e10*
    return 5e-20;
    //return c[0] * std::log(E/c[1])/E*std::exp(-c[2]/E); //fit to ionisation on collision electron-oxygen_1 data obtained via LXcat Source: RBEQ Fit to [Brook et al. J. Phys. B: At. Mol. Phys. 11, 3115 (1978)], [Rothe et al. Phys. Rev. 125, 582 (1962)], [Zipf, Planet Space Sci 33, 1303 (1985)], [Thompson et al. J. Phys. B: At. Mol. Opt. Phys. 28, 1321 (1995)]. Updated: 3 January 2024.
};
void DSMC_MEX_IONIZATION::collide(type_calc3& vel1, type_calc3& vel2, bool& ionize, type_calc& E_new){
    type_calc3 cm = (mass_neus * vel1 + mass_eles * vel2)/sum_mass;

    type_calc3 v_rel = vel1 - vel2;
    type_calc v_rel_mag = v_rel.length();
    #ifdef DEBUG
    type_calc E_ele1 = 0.5 * mass_eles * vel2.length() * vel2.length();
    #endif
    type_calc E_rel = 0.5 * m_reduced * v_rel_mag * v_rel_mag; //in Joules //guess, not physically checked!
    if(E_rel > E_ion_times_10){ // 10* is more phisycally sensible? it means that electron really escapes atom and doesnt feel atom's infuence
        E_new = E_rel - E_ion_times_10;
        ionize = true;
    }
    else{
        ionize = false;
    }

    type_calc cos_ksi = 2*rnd() - 1;
    type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
    type_calc eps = 2*Const::pi * rnd();

    /*rotation*/
    v_rel[0] = v_rel_mag * cos_ksi;
    v_rel[1] = v_rel_mag * sin_ksi*std::cos(eps);
    v_rel[2] = v_rel_mag * sin_ksi*std::sin(eps);

    vel1 = cm + mass_eles/sum_mass*v_rel;
    vel2 = cm - mass_neus/sum_mass*v_rel;

    #ifdef DEBUG
    type_calc E_ele2 = 0.5 * mass_eles * vel2.length() * vel2.length();
    std::cerr << "E_ele1: " << E_ele1 <<  " E_ele2: " << E_ele2 << " E_rel: " << E_rel << " E_ion: " << E_ion_times_10/10 << "\n";
    #endif
    // type_calc E_rel = 0.5 * mass2 * v_rel_mag * v_rel_mag; //in Joules // TRY!!

};


MC_MEX_IONIZATION::MC_MEX_IONIZATION(Species& neutrals, Species& ions, Species& electrons, World& world): neutrals{neutrals}, ions{ions},
 electrons{electrons}, world{world}{ //FOR INTERACTION BETWEEN NEUTRALS AND ELECTRONS, BUT NEUTRALS DONT FEEL ELECTRONS INFLUENCE
                    // Dont erase neutrals suppose that neutrals are so much denser that they do not feel influence of electrons
        if(neutrals.mpw0 < 1e2*electrons.mpw0){
            throw std::invalid_argument("neutrals species should have mpw0 at least 100 times greater than electrons, and density of neutrals' should be much bigger than electrons'");
        }

        if(!neutrals.E_ion || neutrals.E_ion < 0){
            throw std::invalid_argument("neutral species must have proper ionization energy E_ion (in Joules)");
        }
        if(neutrals.mpw0 < ions.mpw0){
            std::cerr << "electrons.mpw0: " << electrons.mpw0 << " ions.mpw0: " << ions.mpw0;
            throw std::invalid_argument("electrons species must have greater or equal mpw0 than ions to prperly calculate number of ions to create during ionisation, best if neutrals.mpw0 = k*ions.mpw0, k in natural numbers");
        }
        E_ion_times_10 = 10*neutrals.E_ion; //It should be not scaled by mpw0 right?
        dv = world.getCellVolume();
        m_reduced = neutrals.mass*electrons.mass/(neutrals.mass + electrons.mass); //TODO SEE BIRD'S APPENDIX A
        mass_neus = neutrals.mass;
        mass_eles = electrons.mass;
        sum_mass = mass_neus + mass_eles;
        num_cells = (world.ni-1) * (world.nj -1) * (world.nk-1);

        c[0] = 1.016e-18;
        c[1] = 9.824e+00;
        c[2] = 6.175e+01;
        c[3] = 666; //not used

		// c[0] = 4.07e-10;
		// c[1] = 0.77;
		// c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		// c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};

type_calc MC_MEX_IONIZATION::evaluateSigma(type_calc v_rel){
    type_calc E = 0.5 * m_reduced * v_rel * v_rel / Const::q_e;
    // type_calc E = 0.5 * electrons.mass * v_rel * v_rel;
    // 1e10*
    return 1e-20;
    //return c[0] * std::log(E/c[1])/E*std::exp(-c[2]/E); //fit to ionisation on collision electron-oxygen_1 data obtained via LXcat Source: RBEQ Fit to [Brook et al. J. Phys. B: At. Mol. Phys. 11, 3115 (1978)], [Rothe et al. Phys. Rev. 125, 582 (1962)], [Zipf, Planet Space Sci 33, 1303 (1985)], [Thompson et al. J. Phys. B: At. Mol. Opt. Phys. 28, 1321 (1995)]. Updated: 3 January 2024.
};

void MC_MEX_IONIZATION::apply(type_calc dt) noexcept{ // test for no pointers
    // std::vector<std::vector<int>> neutrals_in_cell = neutrals.sort_indexes(); // array of array of particle pointers
    // std::vector<std::vector<int>> electrons_in_cell = electrons.sort_indexes(); // array of array of particle pointers
    std::unordered_map<int, std::vector<int>> neutrals_in_cell_index; //will be rare because indexes only for when there are both electrons and neutrons
    std::unordered_map<int, std::vector<int>> electrons_in_cell_index; //rare
    electrons.map_indexes(electrons_in_cell_index); // implement sorting neutrals in map but only for cells that already have electrons
    neutrals.map_indexes(neutrals_in_cell_index, electrons_in_cell_index); //sorts and makes it rarer because only when there are both cpecies in cell
    std::vector<Particle>& parts_neutrals = neutrals.getPartRef();
    std::vector<Particle>& parts_electrons = electrons.getPartRef();

    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = electrons.mpw0;
    int n_collisions = 0;
    int n_ionizations = 0;
    #ifdef DEBUG
    std::cerr << "\n";
    #endif
    for(const std::pair<int, std::vector<int>>& pair: electrons_in_cell_index ){
        int c = pair.first;
        const std::vector<int>& indexes_electrons = pair.second;//electrons_in_cell_index[c];
    
        std::vector<int>& indexes_neutrals = neutrals_in_cell_index[c];
        int np_neus_in_cell = indexes_neutrals.size();
        int np_eles_in_cell = indexes_electrons.size();

        type_calc n_groups_frac = 0.5 * np_neus_in_cell * np_eles_in_cell * neutrals.mpw0 * sigma_v_rel_max * dt / dv; //mpw0 * neutrals.mpw0/mpw0 to scale not equal mpws
        #ifdef DEBUG
        std::cerr << "n_groups_frac: " << n_groups_frac << " c: " << c << " neus: " << np_neus_in_cell << " eles: "<<  np_eles_in_cell << "\n" ;
        #endif
        int n_groups = (int)(n_groups_frac + 0.5);
        if(n_groups > np_neus_in_cell){ //cannot ionize more then it is aviable to ionize
            n_groups = np_neus_in_cell-1;
        }
        //std::cout << "n_groups_frac: " << n_groups_frac << " ";
        for(int g = 0; g < n_groups; g++){
            //std::cout << "np_neus_in_cell: " << np_neus_in_cell << " np_eles_in_cell: " << np_eles_in_cell << "\n";
            // int p_neus = (int)(rnd()*np_neus_in_cell); //sample random neutral particle
            int indexes_neutral_sample_index = rnd(0,np_neus_in_cell);
            int p_neu = indexes_neutrals[indexes_neutral_sample_index]; //sample random neutral particle
            int p_ele = indexes_electrons[rnd(0,np_eles_in_cell)]; //sample random electron particle

            type_calc v_rel = (parts_neutrals[p_neu].vel - parts_electrons[p_ele].vel).length();
            type_calc sigma_v_rel = evaluateSigma(v_rel) * v_rel;
            if(sigma_v_rel > sigma_v_rel_max_temp)
                sigma_v_rel_max_temp = sigma_v_rel;
            type_calc P = sigma_v_rel/sigma_v_rel_max;
            #ifdef DEBUG
            std::cerr << " P: " << P << " sigma_v_rel: " << sigma_v_rel << " sigma_ve_rel_max: " << sigma_v_rel_max << "\n";
            #endif
            bool ionised = false;
                if(P > rnd()){
                n_collisions++;
                type_calc3 vel_new_ele;
                collide(parts_neutrals[p_neu].vel, parts_electrons[p_ele].vel, ionised, vel_new_ele);
                if(ionised){
                    n_ionizations++;
                    int ions_to_create = (int)(electrons.mpw0/ions.mpw0 + 0.5); //assuming electrons.mpw0 >= ions.mpw0
                    for(int i = 0; i < ions_to_create; i++){
                        ions.addParticle(parts_neutrals[p_neu].pos, parts_neutrals[p_neu].vel, ions.mpw0);
                    }

                    // type_calc3 v_rel = parts_neutrals[p_neus]->vel - parts_electrons[p_eles]->vel;
                    // type_calc v_rel_mag = v_rel.length();
                    // type_calc E_rel = 0.5 * neutrals.mass * v_rel_mag * v_rel_mag; //in Joules
                    // type_calc3 lc = world.XtoL(parts_neutrals[p_neus]->pos);
                    // type_calc T = electrons.T.gather(lc);
                    type_calc T = Const::m_e*(vel_new_ele*vel_new_ele)*0.5/Const::k;
                    #ifdef DEBUG
                    std::cerr << "T new electron: " << T << "\n";
                    #endif
                    electrons.addParticle(parts_neutrals[p_neu].pos, vel_new_ele, electrons.mpw0); // change it to random v3th?

                    // Dont erase neutrals suppose that neutrals are so much denser that they do not feel influence of electrons
                }
            }
        }
    }
    if(n_collisions){//change it to acc es only ones to kill
        int np = parts_neutrals.size();
        for(int p = 0; p < np; p++){
            if(parts_neutrals[p].macro_weight == 0){
                parts_neutrals[p] = std::move(parts_neutrals[np-1]);
                p--;
                np--;
            }
        }
        parts_neutrals.erase(parts_neutrals.begin() + np, parts_neutrals.end());
    }
    
    std::cerr << "n_collisions: " << n_collisions << " n_ionizations: " << n_ionizations;
    if(n_collisions){
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
};

void MC_MEX_IONIZATION::collide(const type_calc3& vel1, type_calc3& vel2, bool& ionize, type_calc3& vel_new){
    type_calc3 cm = (mass_neus * vel1 + mass_eles * vel2)/sum_mass;

    type_calc3 v_rel = vel1 - vel2;
    type_calc v_rel_mag = v_rel.length();
    #ifdef DEBUG
    type_calc E_ele1 = 0.5 * mass_eles * vel2.length() * vel2.length();
    #endif
    type_calc E_rel = 0.5 * m_reduced * v_rel_mag * v_rel_mag; //in Joules //guess, not physically checked!
    if(E_rel > E_ion_times_10){ // 10* is more phisycally sensible? it means that electron really escapes atom and doesnt feel atom's infuence
        type_calc E_kin = E_rel - E_ion_times_10;
        ionize = true;
        type_calc r = rnd();
        type_calc E_e1 = E_kin * r;
        type_calc E_new = E_kin * (1-r);
        v_rel_mag = std::sqrt(E_kin*2/m_reduced);
    }
    else{
        ionize = false;
    }

    type_calc cos_ksi = 2*rnd() - 1;
    type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
    type_calc eps = 2*Const::pi * rnd();

    /*rotation*/
    v_rel[0] = v_rel_mag * cos_ksi;
    v_rel[1] = v_rel_mag * sin_ksi*std::cos(eps);
    v_rel[2] = v_rel_mag * sin_ksi*std::sin(eps);

    vel2 = cm - mass_neus/sum_mass*v_rel; //velocity of electrons
    vel_new = cm + mass_eles/sum_mass*v_rel; //velocity of new electrons

    #ifdef DEBUG
    type_calc E_ele2 = 0.5 * mass_eles * vel2.length() * vel2.length();
    std::cerr << "E_ele1: " << E_ele1 <<  " E_ele2: " << E_ele2 << " E_rel: " << E_rel << " E_ion: " << E_ion_times_10/10 << "\n";
    #endif
    // type_calc E_rel = 0.5 * mass2 * v_rel_mag * v_rel_mag; //in Joules // TRY!!

};