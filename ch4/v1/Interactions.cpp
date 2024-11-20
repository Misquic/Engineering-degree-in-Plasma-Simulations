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
                for (int p = 0; p < num_p; p++){
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

    /*delete according number of particles*/
    std::vector<Particle>& particles = neutrals.getPartRef();
    int np = particles.size();
    int deleted = 0;
    type_calc P_sum = 0;
    for(int p = 0; p < np; p ++){
        Particle& part = particles[p];
        type_calc3 lc = world.XtoL(part.pos);
        type_calc dNa = dNi_field.gather(lc);

        int P = dNa*dV/neutrals.mpw0 + rnd();
        P_sum +=P;
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