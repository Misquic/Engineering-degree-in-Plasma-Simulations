#include <algorithm>
#include <unordered_map>
#include <ranges>
#include <fstream>
#include "Interactions.h"
#include "funkc.h"

ChemistryIonize::ChemistryIonize(Species& neutrals, Species& ions, World& world, type_calc rate) noexcept: neutrals{neutrals},
 ions{ions}, world{world}, rate{rate}, dNi_field{world.ni, world.nj, world.nk} {
};
void ChemistryIonize::apply(type_calc dt) noexcept{ //no electrons species //recombination is assumed to only occur on the surface when collision with object
    type_calc dV = world.getCellVolume();
    int created = 0;
    std::vector<std::vector<Particle*>> sorted_neutrals =  neutrals.sortPointers(); // array of array of particle pointers
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
    if(created){
        dmsg("\nions created:" << created);
    }
};


ChemistryIonizeElectrons::ChemistryIonizeElectrons(Species& neutrals, Species& ions, Species& electrons, World& world, type_calc rate) noexcept: neutrals{neutrals},
 ions{ions}, electrons{electrons}, world{world}, rate{rate}, dNi_field{world.ni, world.nj, world.nk} {
};


/*extended version with killing neutrals*/
void ChemistryIonizeElectrons::apply(type_calc dt) noexcept{ //recombination is assumed to only occur on the surface when collision with object
    type_calc dV = world.getCellVolume();
    //dNi_field.clear();
    int created = 0;
    std::vector<std::vector<Particle*>> sorted_neutrals =  neutrals.sortPointers(); // array of array of particle pointers
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
                    type_calc T = electrons.T.gather(lc);
                    ions.addParticle(samped_neutral_ptr->pos, samped_neutral_ptr->vel, ions.mpw0);
                    electrons.addParticle(samped_neutral_ptr->pos, electrons.sampleV3th(T), electrons.mpw0); // change it to random v3th?
                    
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
    
    if(created){
        dmsg("\nions created:" << created);
    }
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
    std::vector<std::vector<Particle*>> parts_in_cell =  species1.sortPointers(); // array of array of particle pointers

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
    if(n_collisions){
        dmsg(" ncollisions: " << n_collisions << "\n");
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
}

void DSMC_MEX::applyTwoSpecies(type_calc dt) noexcept{ //FOR INTERACTION BETWEEN TWO SPECIES
    std::vector<std::vector<Particle*>> parts_in_cell1 =  species1.sortPointers(); // array of array of particle pointers
    std::vector<std::vector<Particle*>> parts_in_cell2 =  species2.sortPointers(); // array of array of particle pointers

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
    if(n_collisions){
        dmsg("ncollisions: " << n_collisions << "\n");
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
DSMC_MEX_Ionization::DSMC_MEX_Ionization(Species& neutrals, Species& ions, Species& electrons, World& world): neutrals{neutrals}, ions{ions},
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


void DSMC_MEX_Ionization::apply(type_calc dt) noexcept{ // test for no pointers
    // std::vector<std::vector<int>> neutrals_in_cell = neutrals.sort_indexes(); // array of array of particle pointers
    // std::vector<std::vector<int>> electrons_in_cell = electrons.sort_indexes(); // array of array of particle pointers
    std::unordered_map<int, std::vector<int>> neutrals_in_cell_index; //will be rare because indexes only for when there are both electrons and neutrons
    std::unordered_map<int, std::vector<int>> electrons_in_cell_index; //rare
    electrons.mapIndexes(electrons_in_cell_index); // implement sorting neutrals in map but only for cells that already have electrons
    neutrals.mapIndexes(neutrals_in_cell_index, electrons_in_cell_index); //sorts and makes it rarer because only when there are both cpecies in cell
    std::vector<Particle>& parts_neutrals = neutrals.getPartRef();
    std::vector<Particle>& parts_electrons = electrons.getPartRef();

    type_calc sigma_v_rel_max_temp = 0;
    type_calc mpw0 = neutrals.mpw0;
    int n_collisions = 0;
    int n_ionizations = 0;

    dmsg("\n");
    
    for(const std::pair<int, std::vector<int>>& pair: electrons_in_cell_index ){
        int c = pair.first;
        const std::vector<int>& indexes_electrons = pair.second;//electrons_in_cell_index[c];
    
        std::vector<int>& indexes_neutrals = neutrals_in_cell_index[c];
        int np_neus_in_cell = indexes_neutrals.size();
        int np_eles_in_cell = indexes_electrons.size();

        type_calc n_groups_frac = 0.5 * np_neus_in_cell * np_eles_in_cell * mpw0 * sigma_v_rel_max * dt / dv;

        dmsg("n_groups_frac: " << n_groups_frac << " c: " << c);

        int n_groups = (int)(n_groups_frac + 0.5);
        if(n_groups > np_neus_in_cell){ //cannot ionize more then it is aviable to ionize
            n_groups = np_neus_in_cell-1;
        }
        for(int g = 0; g < n_groups; g++){
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
            
            dmsg(" P: " << P << " sigma_v_rel: " << sigma_v_rel << " sigma_ve_rel_max: " << sigma_v_rel_max << "\n");

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
                    
                    dmsg("T new electron: " << T << "\n");
                    
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
    if(n_collisions){//change it to access only ones to kill
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
    
    if(n_collisions){
        dmsg("n_collisions: " << n_collisions << " n_ionizations: " << n_ionizations);
        sigma_v_rel_max = sigma_v_rel_max_temp;
    }
};

type_calc DSMC_MEX_Ionization::evaluateSigma(type_calc v_rel){
    type_calc E = 0.5 * m_reduced * v_rel * v_rel / Const::q_e;
    // type_calc E = 0.5 * electrons.mass * v_rel * v_rel;
    // 1e10*
    return 5e-20;
    //return c[0] * std::log(E/c[1])/E*std::exp(-c[2]/E); //fit to ionisation on collision electron-oxygen_1 data obtained via LXcat Source: RBEQ Fit to [Brook et al. J. Phys. B: At. Mol. Phys. 11, 3115 (1978)], [Rothe et al. Phys. Rev. 125, 582 (1962)], [Zipf, Planet Space Sci 33, 1303 (1985)], [Thompson et al. J. Phys. B: At. Mol. Opt. Phys. 28, 1321 (1995)]. Updated: 3 January 2024.
};
void DSMC_MEX_Ionization::collide(type_calc3& vel1, type_calc3& vel2, bool& ionize, type_calc& E_new){
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
    dmsg("E_ele1: " << E_ele1 <<  " E_ele2: " << E_ele2 << " E_rel: " << E_rel << " E_ion: " << E_ion_times_10/10 << "\n");
    #endif
    // type_calc E_rel = 0.5 * mass2 * v_rel_mag * v_rel_mag; //in Joules // TRY!!

};


MC_MEX_Ionization::MC_MEX_Ionization(Species& neutrals, Species& ions, Species& electrons, World& world, std::string coll_data_path, int freq_should_use_map): neutrals{neutrals}, ions{ions},
 electrons{electrons}, world{world}, freq_should_use_map{freq_should_use_map}{ //FOR INTERACTION BETWEEN NEUTRALS AND ELECTRONS, BUT NEUTRALS DONT FEEL ELECTRONS INFLUENCE
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
        E_ion_times_10 = 10*neutrals.E_ion;
        E_ion_J = neutrals.E_ion;
        E_ion_eV = neutrals.E_ion/Const::q_e;
        dv = world.getCellVolume();
        inv_dv = 1/dv;
        dt = world.getDt();
        m_reduced = neutrals.mass*electrons.mass/(neutrals.mass + electrons.mass); //TODO SEE BIRD'S APPENDIX A
        mass_neus = neutrals.mass;
        mass_eles = electrons.mass;
        sum_mass = mass_neus + mass_eles;
        num_cells = (world.ni-1) * (world.nj -1) * (world.nk-1);
        E_precalc_rel_eV = 0.5 * m_reduced/Const::q_e;


    //  a = 1.015e-18, b = 9.793e+00, c = 6.181e+01
        // for evaluateSigmaIon
        c[0] = 1.015e-18;
        c[1] = 9.793e+00;
        c[2] = 6.181e+01;
        c[3] = 666; //not used
        
        B_inc = 10.0; // in eV, for oxygen, argon, but for small energies, couldn't find for bigger source: "A Monte Carlo collision model for the particle-in-cell method: applications to argon and oxygen discharges
                      // V. Vahedi a,b, M. Surendra c
                      // a Department of Electrical Engineering and Computer Science, University of California, Berkeley, CA 94720, USA
                      // b Lawrence Livermore National Laboratory, Livermore, CA 94550, USA
                      // e IBM T.J. Watson Research Center, Yorktown Heights, NY 10598, USA

        // for evaluateSigmaColl
        std::ifstream inputFile;
        inputFile.open(coll_data_path);
        if (!inputFile.is_open()) {
            throw std::invalid_argument("Couldn't open file " + coll_data_path + " for crosssection data of collisions");
        }
        type_calc energy{}, sigma{}, en_max{0}, en_min{std::numeric_limits<type_calc>::max()};
        while(inputFile >> energy >> sigma){
            if(energy > en_max){
                en_max = en_max;
            }else if(energy < en_min){
                en_min = energy;
            }
            data_sigma[energy] = sigma;
            // std::cerr << "input: energy:" << energy << " crosssection: " << sigma << "\n";
        }
        inputFile.close();
        W_sigma_v_rel_max*=std::max(electrons.mpw0, neutrals.mpw0);
		// c[0] = 4.07e-10;
		// c[1] = 0.77;
		// c[2]= 2*Const::k*273.15/m_reduced;	//Bird's reference params at 273.15 K
		// c[3] = std::tgamma(2.5-c[1]); //Gamma(5/2-w)
};

type_calc MC_MEX_Ionization::evaluateSigmaColl(type_calc E){ // expects relative energy in eV
    // type_calc E = 0.5 * electrons.mass * v_rel * v_rel;
    auto upper = data_sigma.lower_bound(E); //interpolation to data from COMMENT: Elastic momentum-transfer cross section | Source: integrated from [WilliamsampAllen J. Phys. B: At. Mol. Opt. Phys. 22, 3529 (1989)] , Above 10 eV: ELSEPA. UPDATED: 2024-01-02 21:39:03
    if(upper == data_sigma.begin()){
        return upper->second;
    }
    if(upper == data_sigma.end()){
        return std::prev(upper)->second;
    }
    auto lower = std::prev(upper);
    type_calc x1 = lower->first;
    type_calc x2 = upper->first;
    type_calc y1 = lower->second;
    type_calc y2 = upper->second;

    type_calc temp = y1 + (E-x1) * (y2-y1) / (x2-x1); 
    return y1 + (E-x1) * (y2-y1) / (x2-x1); 
};
type_calc MC_MEX_Ionization::evaluateSigmaIon(type_calc E){ // expects relative energy in eV
    // type_calc E = E_precalc * v_rel * v_rel;
    // type_calc E = 0.5 * electrons.mass * v_rel * v_rel;
    if(E <= E_ion_eV){
        return 0;
    }
    return c[0] * std::log(E/c[1])/E*std::exp(-c[2]/E); //fit to ionisation on collision electron-oxygen_1 data obtained via LXcat Source: RBEQ Fit to [Brook et al. J. Phys. B: At. Mol. Phys. 11, 3115 (1978)], [Rothe et al. Phys. Rev. 125, 582 (1962)], [Zipf, Planet Space Sci 33, 1303 (1985)], [Thompson et al. J. Phys. B: At. Mol. Opt. Phys. 28, 1321 (1995)]. Updated: 3 January 2024.
};
void MC_MEX_Ionization::apply(type_calc dt) noexcept{
    static bool use_map = true;
    if(world.getTs()%freq_should_use_map == 0){
        std::cerr << "\n";
        use_map = electrons.shouldUseMap();
    }
    if(false){
        dmsg("used map");
        // std::cerr << "\nused map, electron should use map frac: " << electrons.should_use_map();
        // std::cerr << ", size ele: " << electrons.getNumParticles();
        // std::cerr << ", size neus: " << neutrals.getNumParticles();
        // type_calc time_start = world.getWallTime();
        apply_map(dt);   
        // type_calc time_stop = world.getWallTime();
        // std::cerr << ", time: " << time_stop - time_start;


    }else{
        dmsg("used_vector"); 
        // std::cerr << "\nused vector, electron should use map frac: " << electrons.should_use_map();
        // std::cerr << ", size ele: " << electrons.getNumParticles();
        // std::cerr << ", size neus: " << neutrals.getNumParticles();
        // type_calc time_start = world.getWallTime();
        apply_vector_indexes(dt);
        // type_calc time_stop = world.getWallTime();
        // std::cerr << ", time: " << time_stop - time_start;

        //use of pointers is invalid in approach which modifies (adds) particles during appying


    }
}

void MC_MEX_Ionization::apply_vector_indexes(type_calc dt) noexcept{ // using vector when there's a lot of occupied cells using indexes
    std::vector<std::vector<int>> neutrals_in_cell;
    if(neutrals.isSorted()){
        neutrals_in_cell = neutrals.getSortedIndexes(); //already sorted
    }else{
        neutrals_in_cell = neutrals.sortIndexes(); //sort from scratch
    }
    
    std::vector<std::vector<int>> electrons_in_cell;
    if(electrons.isSorted()){
        electrons_in_cell = electrons.getSortedIndexes(); //already sorted
    }else{
        electrons_in_cell = electrons.sortIndexes(); //sort from scratch
    }

    std::vector<Particle>& parts_neutrals = neutrals.getPartRef();
    std::vector<Particle>& parts_electrons = electrons.getPartRef();
    std::vector<Particle>& parts_ions = ions.getPartRef();

    type_calc W_sigma_v_rel_max_temp = 0;
    int n_collisions = 0;
    int n_ionizations = 0;
    int neutrals_to_delete = 0;
    

    dmsg("\n");

    for(int c = 0; c < num_cells; c++){
        std::vector<int>& indexes_electrons = electrons_in_cell[c]; //indexes of electrons in current cell, electrons_in_cell[c]
        int np_eles_in_cell = indexes_electrons.size();
        if(np_eles_in_cell <= 0){
            continue;
        }
        std::vector<int>& indexes_neutrals = neutrals_in_cell[c]; //indexes of neutrals in current cell, neutrals_in_cell[c]
        int np_neus_in_cell = indexes_neutrals.size();
        if(np_neus_in_cell <= 0){
            continue;
        }

        // type_calc neutrals_mpw0 = 0;
        // for(int i: indexes_neutrals){
        //     neutrals_mpw0+=parts_neutrals[i].macro_weight;
        // }
        // neutrals_mpw0/=np_neus_in_cell; //average neutral macroparticle weight in cell //TODO IS IT SENSIBLE? IT BREAKS WHEN NEUTRALS WEIGHT IS LESS Then electrons weight

        //NTC modified
        type_calc n_groups_frac = np_neus_in_cell * np_eles_in_cell * W_sigma_v_rel_max * dt * inv_dv; //mpw0 * neutrals.mpw0/mpw0 to scale not equal mpws, neutrals.mpw0 > electrons.mpw0, 
        
        dmsg("n_groups_frac: " << n_groups_frac << " cell_index: " << c << " np_neus_in_cell: " << np_neus_in_cell << " np_eles_in_cell: "<<  np_eles_in_cell << "\n");
        
        int n_groups = (int)(n_groups_frac + 0.5);
        if(n_groups > np_neus_in_cell){ //cannot ionize more then it is aviable to ionize
            n_groups = np_neus_in_cell-1;
            dmsg("tried to ionize more than neus in cell\n");
        }
        for(int g = 0; g < n_groups; g++){
            int neutral_sample_index = rnd(0,np_neus_in_cell); //index od random neutral in cell
            int p_neu = indexes_neutrals[neutral_sample_index]; //sample random neutral in cell, index from parts_neutrals
            int electron_sample_index = rnd(0,np_eles_in_cell); //index of random electron in cell
            int p_ele = indexes_electrons[electron_sample_index]; //sample random electron in cell, index from parts_electrons
            // type_calc mpw_ele = parts_electrons[p_ele].macro_weight;

            type_calc v_rel = (parts_neutrals[p_neu].vel - parts_electrons[p_ele].vel).length();
            type_calc E_rel = E_precalc_rel_eV * v_rel * v_rel; // in eV
            type_calc sigma_coll = evaluateSigmaColl(E_rel);

            type_calc W_greater = parts_neutrals[p_neu].macro_weight; 
            type_calc W_lesser = parts_electrons[p_ele].macro_weight; 
            greaterLesser(W_greater, W_lesser);
            type_calc W_sigma_v_rel = W_greater* sigma_coll* v_rel;
            if(W_sigma_v_rel > W_sigma_v_rel_max_temp){
                W_sigma_v_rel_max_temp = W_sigma_v_rel;
            }
            type_calc P = W_sigma_v_rel/W_sigma_v_rel_max;
            
            // dmsg(" P: " << P << " sigma_v_rel: " << W_sigma_v_rel << " sigma_ve_rel_max: " << W_sigma_v_rel_max << "\n");
            
            bool ionised = false;
            if(rnd() < P){//variable-weight DSMC
                n_collisions++;
                Particle& electron = parts_electrons[p_ele];
                Particle& neutral = parts_neutrals[p_neu];

                // make created particle could participaate in current sampling 
                if(neutral.macro_weight > electron.macro_weight){ //consider using greater_part, lesser_part and swap, but check if it swaps only pointers or particles
                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //split neutral, differential weight in place of previous, weight of ele in new place
                    neutral.macro_weight-= electron.macro_weight; //serves as deleting of neutral
                    //apply collision
                    type_calc3 vel_new_ele{};
                    type_calc3 vel_splitted_neu = neutral.vel;
                    collide(vel_splitted_neu, electron.vel, ionised, vel_new_ele, sigma_coll);
                    if(ionised){
                        n_ionizations++;
                        parts_ions.emplace_back(neutral.pos, vel_splitted_neu, W_lesser); //without velocity rewind
                        parts_electrons.emplace_back(neutral.pos, vel_new_ele, W_lesser); //without velocity rewind
                        
                    }else{
                        //rest of split// If I use event splitting, I will need to split every time so it will be more streamlined 
                        parts_neutrals.emplace_back(neutral.pos, neutral.vel, electron.macro_weight); //new neutral of weight of electron
                        indexes_neutrals.push_back(parts_neutrals.size()-1); //add index of new neutral to indexes in cell
                        np_neus_in_cell++;

                    }
                    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                }else if(neutral.macro_weight < electron.macro_weight){
                    //split electron, differential weight in new place, weight of neu in current place
                    electron.macro_weight=neutral.macro_weight; //splited electron

                    parts_electrons.emplace_back(electron.pos, electron.vel, electron.macro_weight - neutral.macro_weight); // differential electron
                    indexes_electrons.push_back(parts_electrons.size()-1); //add index of differential electron to indexes in cell

                    type_calc3 vel_new_ele{};
                    collide(neutral.vel, electron.vel, ionised, vel_new_ele, sigma_coll);
                    if(ionised){
                        n_ionizations++;
                        parts_ions.emplace_back(neutral.pos, neutral.vel, W_lesser); //without velocity rewind
                        parts_electrons.emplace_back(neutral.pos, vel_new_ele, W_lesser); //without velocity rewind
                        //delete neutral
                        np_neus_in_cell--;
                        indexes_neutrals[neutral_sample_index] = indexes_neutrals[np_neus_in_cell]; // use of previous decrement, if used wwithout previous decrement then should be np_neus_in_cell-1
                        parts_neutrals[p_neu] = std::move(parts_neutrals.back()); 
                        neutrals_to_delete++;
                        //TODO mark neutrals as nonsorted

                    }

                }else{//equal weights
                                    // type_calc3 vel_new_ele;
                // collide(parts_neutrals[p_neu].vel, parts_electrons[p_ele].vel, ionised, vel_new_ele, sigma_coll);
                // if(ionised){
                //     n_ionizations++;
                //     ions.addParticle(parts_neutrals[p_neu].pos, parts_neutrals[p_neu].vel, mpw_ele);
                //     electrons.addParticle(parts_neutrals[p_neu].pos, vel_new_ele, mpw_ele);
                }


                // type_calc3 vel_new_ele;
                // collide(parts_neutrals[p_neu].vel, parts_electrons[p_ele].vel, ionised, vel_new_ele, sigma_coll);
                // if(ionised){
                //     n_ionizations++;
                //     ions.addParticle(parts_neutrals[p_neu].pos, parts_neutrals[p_neu].vel, mpw_ele);
                //     electrons.addParticle(parts_neutrals[p_neu].pos, vel_new_ele, mpw_ele);

                //     // Dont erase neutrals suppose that neutrals are so much denser that they do not feel influence of electrons
                // }
            }
        }
    }


    if(n_collisions){
        ions.setSorted(false);
        electrons.setSorted(false);
        neutrals.setSorted(false);
        dmsg("n_collisions: " << n_collisions << " n_ionizations: " << n_ionizations << " ");
        W_sigma_v_rel_max = W_sigma_v_rel_max_temp;
        if(neutrals_to_delete){
            parts_neutrals.erase(parts_neutrals.begin() + parts_neutrals.size() - neutrals_to_delete, parts_neutrals.end());
        }
    }
    dmsg("interaction end\n");
}
void MC_MEX_Ionization::apply_map(type_calc dt) noexcept{ // using map when there's NOT a lot of occupied cells
    // // std::vector<std::vector<int>> neutrals_in_cell = neutrals.sort_indexes(); // array of array of particle pointers
    // // std::vector<std::vector<int>> electrons_in_cell = electrons.sort_indexes(); // array of array of particle pointers
    // std::unordered_map<int, std::vector<int>> neutrals_in_cell_index; //will be rare because indexes only for when there are both electrons and neutrons
    // std::unordered_map<int, std::vector<int>> electrons_in_cell_index; //rare
    // electrons.mapIndexes(electrons_in_cell_index); // implement sorting neutrals in map but only for cells that already have electrons
    // neutrals.mapIndexes(neutrals_in_cell_index, electrons_in_cell_index); //sorts and makes it rarer because only when there are both cpecies in cell
    // std::vector<Particle>& parts_neutrals = neutrals.getPartRef();
    // std::vector<Particle>& parts_electrons = electrons.getPartRef();

    // type_calc sigma_v_rel_max_temp = 0;
    // int n_collisions = 0;
    // int n_ionizations = 0;

    // dmsg("\n");
    
    // for(const std::pair<int, std::vector<int>>& pair: electrons_in_cell_index ){
    //     int c = pair.first;
    //     const std::vector<int>& indexes_electrons = pair.second;//electrons_in_cell_index[c];
    
    //     std::vector<int>& indexes_neutrals = neutrals_in_cell_index[c];
    //     int np_neus_in_cell = indexes_neutrals.size();
    //     if(np_neus_in_cell <= 0){
    //         continue;
    //     }
    //     int np_eles_in_cell = indexes_electrons.size();

    //     type_calc neutrals_mpw0 = 0;
    //     for(int i: indexes_neutrals){
    //         neutrals_mpw0+=parts_neutrals[i].macro_weight;
    //     }
    //     neutrals_mpw0/=np_neus_in_cell; //average neutral macroparticle weight in cell //TODO IS IT SENSIBLE? IT BREAKS WHEN NEUTRALS WEIGHT IS LESS Then electrons weight

    //     type_calc n_groups_frac = 0.5 * np_neus_in_cell * np_eles_in_cell * neutrals_mpw0 * sigma_v_rel_max * dt / dv; //mpw0 * neutrals.mpw0/mpw0 to scale not equal mpws, neutrals.mpw0 > electrons.mpw0, 
        
    //     dmsg("n_groups_frac: " << n_groups_frac << " c: " << c << " neus: " << np_neus_in_cell << " eles: "<<  np_eles_in_cell << "\n");
        
    //     int n_groups = (int)(n_groups_frac + 0.5);
    //     if(n_groups > np_neus_in_cell){ //cannot ionize more then it is aviable to ionize
    //         n_groups = np_neus_in_cell-1;
    //     }
    //     for(int g = 0; g < n_groups; g++){
    //         int indexes_neutral_sample_index = rnd(0,np_neus_in_cell);
    //         int p_neu = indexes_neutrals[indexes_neutral_sample_index]; //sample random neutral particle
    //         int p_ele = indexes_electrons[rnd(0,np_eles_in_cell)]; //sample random electron particle
    //         type_calc mpw_ele = parts_electrons[p_ele].macro_weight;

    //         type_calc v_rel = (parts_neutrals[p_neu].vel - parts_electrons[p_ele].vel).length();
    //         type_calc E_rel = E_precalc_rel_eV * v_rel * v_rel; // in eV
    //         type_calc sigma_coll = evaluateSigmaColl(E_rel);
    //         type_calc sigma_v_rel = sigma_coll* v_rel;
    //         if(sigma_v_rel > sigma_v_rel_max_temp)
    //             sigma_v_rel_max_temp = sigma_v_rel;
    //         type_calc P = sigma_v_rel/sigma_v_rel_max;
            
    //         //dmsg(" P: " << P << " sigma_v_rel: " << sigma_v_rel << " sigma_ve_rel_max: " << sigma_v_rel_max << "\n");
            
    //         bool ionised = false;
    //         if(P > rnd()){
    //             n_collisions++;
    //             type_calc3 vel_new_ele;
    //             collide(parts_neutrals[p_neu].vel, parts_electrons[p_ele].vel, ionised, vel_new_ele, sigma_coll);
    //             if(ionised){
    //                 n_ionizations++;
    //                 ions.addParticle(parts_neutrals[p_neu].pos, parts_neutrals[p_neu].vel, mpw_ele);
    //                 electrons.addParticle(parts_neutrals[p_neu].pos, vel_new_ele, mpw_ele);

    //                 // Dont erase neutrals suppose that neutrals are so much denser that they do not feel influence of electrons
    //             }
    //         }
    //     }
    // }
    
    // if(n_collisions){
    //     ions.setSorted(false);
    //     electrons.setSorted(false);
    //      dmsg("n_collisions: " << n_collisions << " n_ionizations: " << n_ionizations << " ");
    //     sigma_v_rel_max = sigma_v_rel_max_temp;
    // }
    // dmsg("interaction end\n");
};

void MC_MEX_Ionization::newVelocityElecton(type_calc3& vel, type_calc E, type_calc3& vel_inc_unit){ // in LAB
        type_calc cos_ksi = (2 + E - 2*std::pow( 1+E, rnd())) / E;
        type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
        type_calc phi =  2*Const::pi * rnd();
        //type_calc inv_sin_theta = 1/(std::sqrt(1 - std::pow(vel_inc_unit*type_calc3{1,0,0}, 2)));
        type_calc v_mag = std::sqrt(E*two_qe_to_me);
        #ifdef DEBUG
        if(std::isnan(v_mag)){
            std::cerr<< "nan new vel ele v_mag\n";
            v_mag = 0;
        }
        #endif

        static type_calc3 i(1,0,0);
        type_calc3 i_cross_v_inc = i.cross(vel_inc_unit);
        vel = cos_ksi*vel_inc_unit + i_cross_v_inc*sin_ksi*std::sin(phi) + vel_inc_unit.cross(i_cross_v_inc)*sin_ksi*std::cos(phi); 
        vel *= v_mag;
   
        // vel = vel_inc_unit*v_mag;
        // vel.rotate(vel_inc_unit.cross({1,0,0}).unit(), sin_ksi, cos_ksi);
        // vel.rotate(vel_inc_unit, std::sin(phi), std::cos(phi));


        // vel = vel_inc_unit*cos_ksi + v_inc_cross_i * (sin_ksi * std::sin(phi)*inv_sin_theta) - vel_inc_unit.cross(v_inc_cross_i) * sin_ksi*cos_ksi*inv_sin_theta;
        // vel *= v_mag;
};
void MC_MEX_Ionization::newVelocityElecton(type_calc3& vel, type_calc E){ // in CM
    type_calc v_scattered_mag = std::sqrt(E*two_qe_to_me);

    type_calc cos_ksi = 2*rnd() - 1;
    type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
    type_calc eps = 2*Const::pi * rnd();

    /*rotation*/
    vel[0] =v_scattered_mag * cos_ksi;
    vel[1] =v_scattered_mag * sin_ksi*std::cos(eps);
    vel[2] =v_scattered_mag * sin_ksi*std::sin(eps);
};

#define IONIZE_1 // IONIZE_1 - collision with ionisation in LAB || IONIZE_2 - in CM
void MC_MEX_Ionization::collide(type_calc3& vel_neu, type_calc3& vel_ele, bool& ionize, type_calc3& vel_new, type_calc sigma_coll){
    type_calc3 v_rel = vel_neu - vel_ele;
    type_calc v_rel_mag = v_rel.length();
    #ifdef DEBUG
    type_calc E_ele1 = 0.5 * mass_eles * vel_ele.length() * vel_ele.length();
    type_calc sigma = evaluateSigmaIon(v_rel_mag);
    #endif
    type_calc E_rel_J = 0.5 * m_reduced * v_rel_mag * v_rel_mag; //in Joules
    // type_calc P = 1-std::exp(-N_atoms_macro*evaluateSigmaIon(v_rel_mag)*v_rel_mag*dt*inv_dv);
    type_calc P = evaluateSigmaIon(E_rel_J/Const::q_e)/sigma_coll; 
    if(rnd() <= P){ 
        ionize = true; 
        type_calc vel_ele_mag = vel_ele.length();

        #ifdef IONIZE_1 // collision with ionisation in LAB
        type_calc E_incident_el = vel_ele_mag*vel_ele_mag*E_precalc_ele_eV; //precalc used for faster math
        if(E_incident_el < E_ion_eV){
            ionize = false;
            std::cerr << "E_incident_el - E_ion_eV < 0: " << E_incident_el - E_ion_eV << "\n";
            return;
        }
        type_calc E_ejected = 10.0*std::tan(rnd()*std::atan((E_incident_el - E_ion_eV)/(2*B_inc)));
        type_calc E_scattered = E_incident_el - E_ion_eV - E_ejected;
        if(E_scattered < 0) E_scattered = 0.000001;
        //scattered electron:
        type_calc3 vel_inc_unit = vel_ele.unit();
        newVelocityElecton(vel_ele, E_scattered, vel_inc_unit);
        newVelocityElecton(vel_new, E_ejected, vel_inc_unit);

        #else
        #ifdef IONIZE_2 //in CM
        type_calc3 cm = (mass_neus * vel_neu + mass_eles * vel_ele)/sum_mass;

        type_calc E_kin = E_rel_J - E_ion_J;

        type_calc E_scattered = E_kin * rnd();
        newVelocityElecton(v_rel, E_scattered);
        
        vel_ele = cm - mass_neus/sum_mass*v_rel;

        type_calc E_ejected = E_kin - E_scattered;
        newVelocityElecton(v_rel, E_ejected);

        vel_new = cm + mass_eles/sum_mass*v_rel;
        #endif
        #endif

    }
    else{
        ionize = false; // normal collision ion doesn't feels impact, in CM
        type_calc3 cm = (mass_neus * vel_neu + mass_eles * vel_ele)/sum_mass;

        type_calc cos_ksi = 2*rnd() - 1;
        type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
        type_calc eps = 2*Const::pi * rnd();

        /*rotation*/
        v_rel[0] = v_rel_mag * cos_ksi;
        v_rel[1] = v_rel_mag * sin_ksi*std::cos(eps);
        v_rel[2] = v_rel_mag * sin_ksi*std::sin(eps);

        vel_ele = cm - mass_neus/sum_mass*v_rel; //velocity of old electron
    }
    ////////////////////////////////////////////////
    // if(rnd() <= P){
    // // if(E_rel > E_ion_times_10/10.0){ // 10* is more phisycally sensible? it means that electron really escapes atom and doesnt feel atom's infuence
    //     type_calc E_kin = E_rel - E_ion;
    //     ionize = true;
    //     type_calc r = rnd();
    //     type_calc E_e1 = E_kin * r;
    //     type_calc E_new = E_kin * (1-r);
    //     v_rel_mag = std::sqrt(E_kin*2/m_reduced);
    // }
    // else{
    //     ionize = false;
    // }

    // type_calc cos_ksi = 2*rnd() - 1;
    // type_calc sin_ksi = std::sqrt(1-cos_ksi*cos_ksi);
    // type_calc eps = 2*Const::pi * rnd();

    // /*rotation*/
    // v_rel[0] = v_rel_mag * cos_ksi;
    // v_rel[1] = v_rel_mag * sin_ksi*std::cos(eps);
    // v_rel[2] = v_rel_mag * sin_ksi*std::sin(eps);

    // vel2 = cm - mass_neus/sum_mass*v_rel; //velocity of old electron
    // vel_new = cm + mass_eles/sum_mass*v_rel; //velocity of new electrons

    // #ifdef DEBUG
    // type_calc E_ele2 = 0.5 * mass_eles * vel2.length() * vel2.length();
    // std::cerr << "E_ele1: " << E_ele1 <<  " E_ele2: " << E_ele2 << " E_rel: " << E_rel << " E_ion: " << E_ion_times_10/10 << "\n";
    // #endif
    // type_calc E_rel = 0.5 * mass2 * v_rel_mag * v_rel_mag; //in Joules // TRY!!

};