#include "Species.h"
#include "Rnd.h"


/*Particle constructors*/
Particle::Particle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight) noexcept: pos{x, y, z}, vel{u, v, w}, macro_weight{macro_weight} {
};
Particle::Particle(type_calc3 pos, type_calc3 vel, type_calc macro_weight) noexcept: pos{pos}, vel{vel}, macro_weight{macro_weight} {
};

/*Species constructors*/
Species::Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0): name{name}, mass{mass},
 charge{charge}, world{world}, mpw0{mpw0}, den{world.nn}, den_avg{world.nn}, n_sum{world.nn}, nv_sum{world.nn}, nuu_sum{world.nn},
 nvv_sum{world.nn}, nww_sum{world.nn}, T{world.nn}, vel{world.nn}, macro_part_count{world.ni -1, world.nj -1, world.nk-1} {
    sorted_particles_pointers.reserve((world.ni-1)*(world.nj-1)*(world.nk-1));
};
Species::Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0, type_calc E_ion): name{name}, mass{mass},
 charge{charge}, world{world}, mpw0{mpw0}, E_ion{E_ion}, den{world.nn}, den_avg{world.nn}, n_sum{world.nn}, nv_sum{world.nn}, nuu_sum{world.nn},
 nvv_sum{world.nn}, nww_sum{world.nn}, T{world.nn}, vel{world.nn}, macro_part_count{world.ni -1, world.nj -1, world.nk-1} {
    sorted_particles_pointers.reserve((world.ni-1)*(world.nj-1)*(world.nk-1));
};


/*Species methods*/
size_t Species::getNumParticles() const{
    return particles.size();
};
void Species::advanceSputtering(Species& neutrals, Species& spherium, type_calc dt){//TODO change it so species take into account if its neutral or not? -> subclasses
    type_calc3 x0 = world.getX0();
    type_calc3 xm = world.getXm();

    type_calc3 lc{};
    type_calc3 ef_part{};
    
    size_t np = particles.size();

    for(size_t p = 0; p < np; p++){
        Particle& part = particles[p];
        type_calc3 lc = world.XtoL(part.pos);

        type_calc3 ef_part = world.ef.gather(lc);

        part.vel += ef_part*(dt*charge/mass);

        type_calc t_reminding = 1;
        int n_bounces = 0;
        while(t_reminding > 0){
            if(++n_bounces > 20) {
                particles[p] = std::move(particles[np-1]);
                np--;
                p--;
                break;
            }

            type_calc3 pos_old = part.pos;
            part.pos += part.vel*t_reminding*dt;

            int in_object = world.inObject(part.pos);

            if(!world.inBounds(part.pos)){
                particles[p] = std::move(particles[np-1]);
                np--;
                p--;
                break;
            }
            else if(in_object){
                type_calc tp;
                type_calc3 n; //normal vector to the surface at the point of intersection
                world.lineIntersect(pos_old, part.pos, in_object, tp, part.pos, n); // passing tp, pos, n so it overwrites them to correct values, in stead of returning tuple
                type_calc v_mag_pre_impact = part.vel.length();
                if(!charge){ //neutrals
                    part.vel = sampleReflectedVelocity(part.pos, part.vel.length(), n); 
                }
                else{ //ions
                    type_calc mpw_ratio = this->mpw0/neutrals.mpw0;

                    /*inject neutrals*/
                    int mp_create = (int)(mpw_ratio + rnd()); //number of macroparticles to create
                    for(int i = 0; i < mp_create; i ++){
                        type_calc3 vel = sampleReflectedVelocity(part.pos, v_mag_pre_impact, n);
                        neutrals.addParticle(part.pos, vel);
                    }
                    /*emit sputtered material using simple yield model*/
                    type_calc sput_yield = (v_mag_pre_impact>5e3)?0.1:0;
                    type_calc sput_macroweight_ratio = sput_yield*this->mpw0/spherium.mpw0;
                    int sput_mp_create = (int)(sput_macroweight_ratio + rnd());//number of sputtered macroparticles to create
                    for(int i = 0; i < sput_mp_create; i++){
                        type_calc3 vel = sampleReflectedVelocity(part.pos, v_mag_pre_impact, n);
                        spherium.addParticle(part.pos, vel);
                    }
                    particles[p] = std::move(particles[np-1]);
                    np--;
                    p--;
                    break;
                }

                t_reminding *= (1-tp);
                continue;
            }
            
            t_reminding = 0;
        }
    }
    particles.erase(particles.begin() + np, particles.end());
};
//
void Species::advanceNoSputtering(Species& neutrals, type_calc dt){//TODO change it so species take into account if its neutral or not? -> subclasses
    type_calc3 x0 = world.getX0();
    type_calc3 xm = world.getXm();

    type_calc3 lc{};
    type_calc3 ef_part{};
    
    size_t np = particles.size();

    for(size_t p = 0; p < np; p++){
        Particle& part = particles[p];
        type_calc3 lc = world.XtoL(part.pos);

        type_calc3 ef_part = world.ef.gather(lc);

        part.vel += ef_part*(dt*charge/mass);

        type_calc t_reminding = 1;
        int n_bounces = 0;
        while(t_reminding > 0){
            if(++n_bounces > 20) {
                particles[p] = std::move(particles[np-1]);
                np--;
                p--;
                break;
            }

            type_calc3 pos_old = part.pos;
            part.pos += part.vel*t_reminding*dt;

            int in_object = world.inObject(part.pos);

            if(!world.inBounds(part.pos)){
                particles[p] = std::move(particles[np-1]);
                np--;
                p--;
                break;
            }
            else if(in_object){
                type_calc tp;
                type_calc3 n; //normal vector to the surface at the point of intersection
                world.lineIntersect(pos_old, part.pos, in_object, tp, part.pos, n); // passing tp, pos, n so it overwrites them to correct values, in stead of returning tuple
                type_calc v_mag_pre_impact = part.vel.length();
                if(!charge){ //neutrals
                    part.vel = sampleReflectedVelocity(part.pos, part.vel.length(), n); 
                }
                else{ //ions
                    type_calc mpw_ratio = this->mpw0/neutrals.mpw0;

                    /*inject neutrals*/
                    int mp_create = (int)(mpw_ratio + rnd()); //number of macroparticles to create
                    for(int i = 0; i < mp_create; i ++){
                        type_calc3 vel = sampleReflectedVelocity(part.pos, v_mag_pre_impact, n);
                        neutrals.addParticle(part.pos, vel);
                    }                    
                    particles[p] = std::move(particles[np-1]);
                    np--;
                    p--;
                    break;
                }

                t_reminding *= (1-tp);
                continue;
            }
            
            t_reminding = 0;
        }
    }
    particles.erase(particles.begin() + np, particles.end());
};
void Species::advanceElectrons(type_calc dt){
    type_calc3 x0 = world.getX0();
    type_calc3 xm = world.getXm();

    type_calc3 lc{};
    type_calc3 ef_part{};
    
    size_t np = particles.size();

    for(size_t p = 0; p < np; p++){
        Particle& part = particles[p];
        type_calc3 lc = world.XtoL(part.pos);

        type_calc3 ef_part = world.ef.gather(lc);

        part.vel += ef_part*(dt*charge/mass);
        part.pos += part.vel*dt;

        int in_object = world.inObject(part.pos);

        if(!world.inBounds(part.pos)){
            particles[p] = std::move(particles[np-1]);
            np--;
            p--;
            continue;
        }
        else if(in_object){                    
            particles[p] = std::move(particles[np-1]);
            np--;
            p--;
            continue;
        }
    }
    particles.erase(particles.begin() + np, particles.end());
};
void Species::computeNumberDensity(){ //to ommit doind this thice implement using n_sum

    den.clear();
    for(Particle& part : particles)
    {
        type_calc3 lc = world.XtoL(part.pos);
        den.scatter(lc, part.macro_weight);
    }

    den/=world.node_vol;
};
void Species::addParticle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight){
    addParticle({x, y, z}, {u, v, w}, macro_weight);
};
void Species::addParticle(type_calc3 pos, type_calc3 vel, type_calc macro_weight){
    if (!world.inBounds(pos)) return;
    if(world.inObject(pos)) return;

    type_calc3 lc = world.XtoL(pos);

    type_calc3 ef_part = world.ef.gather(lc);

    vel -= charge/mass*ef_part*(0.5 * world.getDt());

    particles.emplace_back(pos, vel, macro_weight);
};
void Species::addParticle(type_calc3 pos, type_calc3 vel){
    addParticle(pos, vel, mpw0);
};

void Species::loadParticleBox(type_calc3 x1, type_calc3 x2, type_calc num_den_, type_calc num_macro){
    throw std::runtime_error("Species::loadParticleBox - use of depracated method");
    
    /*x1 - starting vertex of box
      x2 - opposite vertex of box
      num den - 
      num_macro_particles - number of macro particles*/
    type_calc box_vol = (x2[0]-x1[0])*(x2[1]-x1[1])*(x2[2]-x1[2]);
    type_calc num_micro = num_den_ * box_vol;
    type_calc macro_weight = num_micro/num_macro; 
    // this->mpw0 = macro_weight;

    std::cout << "Loading particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weigh of macroparticle: " << macro_weight << "\n";

    particles.reserve(num_macro);

    type_calc3 pos;
    type_calc3 vel;
    for(int i = 0; i < num_macro; i++){
        pos = {rnd(x1[0], x2[0]), rnd(x1[1], x2[1]), rnd(x1[2], x2[2])};
        vel = {0, 0, 0}; //stationary particle;
        // if(i > 50 && i < 55){
        //     std::cout << pos << "\n";
        // }

        addParticle(pos, vel, macro_weight);
    }
};
void Species::loadParticleBoxQS(type_calc3 xc, type_calc3 sides, type_calc num_den_, int3 num_macro){ //QS - quiet start
    /*x1 - starting vertex of box
      x2 - opposite vertex of box
      num den - 
      num_macro - number of macro particles per axis*/
    type_calc box_vol = sides[0]*sides[1]*sides[2];
    int num_macro_tot = (num_macro[0]-1)*(num_macro[1]-1)*(num_macro[2]-1);
    type_calc num_micro = num_den_ * box_vol;
    type_calc macro_weight = num_micro/num_macro_tot; 

    std::cout << "Quiet Start, loading (uniform spacing) particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << "\nwith number of macroparticles: " << num_macro_tot << " with grid of macroparticles: " << num_macro << " and mass of macroparticle: " << macro_weight << "\n\n";

    particles.reserve(num_macro_tot);

    //particle grid spacing
    type_calc di = (sides[0])/(num_macro[0]-1);
    type_calc dj = (sides[1])/(num_macro[1]-1);
    type_calc dk = (sides[2])/(num_macro[2]-1);

    type_calc3 x1 = xc-sides/2;
    type_calc3 x2 = xc+sides/2;

    type_calc3 pos{};
    type_calc3 vel{};
    type_calc w{}; //number weight for mass of macroparticle to counterpart not equal cell volumes

    for(int i = 0; i < num_macro[0]; i++){
        for(int j = 0; j < num_macro[1]; j++){
            for(int k = 0; k < num_macro[2]; k++){

                pos = {x1[0] + i*di,
                       x1[1] + j*dj, 
                       x1[2] + k*dk};
                vel = {0, 0, 0}; //stationary particle;

                //shifting of particles that are on max faces of grid back to domain
                if(pos[0] >= x2[0]) pos[0] -= 5e-4*di;
                if(pos[1] >= x2[1]) pos[1] -= 5e-4*dj;
                if(pos[2] >= x2[2]) pos[2] -= 5e-4*dk;

                w = 1;
                if(i == 0 || i == num_macro[0]-1) w*=0.5;
                if(j == 0 || j == num_macro[1]-1) w*=0.5;
                if(k == 0 || k == num_macro[2]-1) w*=0.5; 

                addParticle(pos, vel, macro_weight*w);
            }
        }
    }
};

void Species::loadParticleBoxQS(type_calc3 xc, type_calc3 sides, type_calc num_den_){ //QS - quiet start
    /*x1 - starting vertex of box
      x2 - opposite vertex of box
      num den - 
      num_macro - number of macro particles per axis*/
    type_calc box_vol = sides[0]*sides[1]*sides[2];
    // int num_macro_tot = (num_macro[0]-1)*(num_macro[1]-1)*(num_macro[2]-1);
    type_calc num_micro = num_den_ * box_vol;
    int num_macro_tot = int(num_micro/mpw0);
    int grid_spacing = std::cbrt(num_macro_tot);
    int3 grid(grid_spacing);
    num_macro_tot = std::pow(grid_spacing,3);

    std::cout << "Quiet Start, loading (uniform spacing) particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << "\nwith number of macroparticles: " << num_macro_tot << " with grid of macroparticles: " << grid << " and mass of macroparticle: " << mpw0<< "\n\n";
    std::cout << std::flush;
    particles.reserve(num_macro_tot);

    //particle grid spacing
    type_calc di = (sides[0])/(grid[0]-1);
    type_calc dj = (sides[1])/(grid[1]-1);
    type_calc dk = (sides[2])/(grid[2]-1);

    type_calc3 x1 = xc-sides/2;
    type_calc3 x2 = xc+sides/2;

    type_calc3 pos{};
    type_calc3 vel{};
    type_calc w{}; //number weight for mass of macroparticle to counterpart not equal cell volumes

    for(int i = 0; i < grid[0]; i++){
        for(int j = 0; j < grid[1]; j++){
            for(int k = 0; k < grid[2]; k++){

                pos = {x1[0] + i*di,
                       x1[1] + j*dj, 
                       x1[2] + k*dk};
                vel = {0, 0, 0}; //stationary particle;

                //shifting of particles that are on max faces of grid back to domain
                if(pos[0] >= x2[0]) pos[0] -= 5e-4*di;
                if(pos[1] >= x2[1]) pos[1] -= 5e-4*dj;
                if(pos[2] >= x2[2]) pos[2] -= 5e-4*dk;

                addParticle(pos, vel, mpw0);
            }
        }
    }
};

void Species::loadParticleBoxThermal(type_calc3 x0, type_calc3 sides, type_calc num_den_, type_calc T){
    type_calc box_vol = sides[0]*sides[1]*sides[2];
    type_calc num_micro = num_den_ * box_vol;
    size_t num_macro = (size_t)(num_micro/mpw0); 
    // this->mpw0 = macro_weight;

    std::cout << "Loading particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weigh of macroparticle: " << this->mpw0 << std::endl;
    if(num_macro < 1){
        throw std::invalid_argument("number of macroparticles less than 1, change initial values\n");
    }

    this->particles.reserve(num_macro);

    type_calc3 pos;
    type_calc3 vel;
    type_calc3 xmin = x0-sides/2;
    type_calc3 xmax = x0+sides/2;
    int every_1_pm = num_macro/100;
    if(every_1_pm < 1){
        every_1_pm = 1;
    }else{
        while(every_1_pm > 1e6){
            every_1_pm /= 10;
        }
    }
    for(size_t i = 0; i < num_macro; i++){
        if(i%every_1_pm==0){
            std::cout << "\r                                         \rloaded: " << (type_calc)i/num_macro*100.0 << "%" << std::flush;
        }
        pos = {rnd(xmin[0], xmax[0]), rnd(xmin[1], xmax[1]), rnd(xmin[2], xmax[2])};
        vel = sampleV3th(T);
        addParticle(pos, vel, mpw0);
    }
    std::cout << "\rLoaded number of macroparticles: " << particles.size()<< "\n";

};

void Species::loadParticleSphereThermal(type_calc3 x0, type_calc r, type_calc num_den_, type_calc T){
    type_calc box_vol = 4/3*Const::pi*r*r*r;
    type_calc num_micro = num_den_ * box_vol;
    size_t num_macro = (size_t)(num_micro/mpw0); 
    // this->mpw0 = macro_weight;

    std::cout << "Loading particle " << name << " in a sphere of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weigh of macroparticle: " << this->mpw0 << std::endl;
    if(num_macro < 1){
        throw std::invalid_argument("number of macroparticles less than 1, change initial values\n");
    }

    this->particles.reserve(num_macro);

    type_calc3 pos;
    type_calc3 vel;

    int every_1_pm = num_macro/100;
    if(every_1_pm < 1){
        every_1_pm = 1;
    }else{
        while(every_1_pm > 1e6){
            every_1_pm /= 10;
        }
    }
    for(size_t i = 0; i < num_macro; i++){
        if(i%every_1_pm==0){
            std::cout << "\r                                         \rloaded: " << (type_calc)i/num_macro*100.0 << "%" << std::flush;
        }
        type_calc r_rnd = r*std::cbrt(rnd());
        type_calc phi = 2*Const::pi*rnd();
        type_calc sin_phi = std::sin(phi);
        type_calc cos_phi = std::cos(phi);
        type_calc cos_theta = 2*rnd() - 1;
        type_calc sin_theta = std::sqrt(1-cos_theta*cos_theta);
        pos = x0 + type_calc3{r_rnd*sin_theta*cos_phi, r_rnd*sin_theta*sin_phi, r_rnd*cos_theta};
        vel = sampleV3th(T);
        addParticle(pos, vel, mpw0);
    }
    std::cout << "\rLoaded number of macroparticles: " << particles.size()<< "\n";

};

void Species::loadLangumir(type_calc lambda, type_calc A, type_calc num_den_, type_calc T){
    type_calc box_vol = world.getL()[0]*world.getL()[1]*world.getL()[2];
    type_calc num_micro = num_den_ * box_vol;
    size_t num_macro = (size_t)(num_micro/mpw0);
    type_calc k = 2*Const::pi/lambda;
    // this->mpw0 = macro_weight;

    std::cout << "Loading Lagumir " << name << " in a whole domain of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weigh of macroparticle: " << this->mpw0 << std::endl;
    if(num_macro < 1){
        throw std::invalid_argument("number of macroparticles less than 1, change initial values\n");
    }

    this->particles.reserve(num_macro);

    type_calc3 pos;
    type_calc3 vel;

    type_calc3 xmin = world.getXc()-world.getL()/2;
    type_calc3 xmax = world.getXc()+world.getL()/2;

    int every_1_pm = num_macro/100;
    if(every_1_pm < 1){
        every_1_pm = 1;
    }else{
        while(every_1_pm > 1e6){
            every_1_pm /= 10;
        }
    }
    for(size_t i = 0; i < num_macro; i++){
        if(i%every_1_pm==0){
            std::cout << "\r                                         \rloaded: " << (type_calc)i/num_macro*100.0 << "%" << std::flush;
        }
        type_calc x = rnd(xmin[0], xmax[0]);
        pos = {x+ A*std::sin(k*x), rnd(xmin[1], xmax[1]), rnd(xmin[2], xmax[2])};
        vel = sampleV3th(T);
        addParticle(pos, vel, mpw0);
    }
    std::cout << "\rLoaded number of macroparticles: " << particles.size()<< "\n";

};
void Species::move(type_calc3 dx, Rectangle& rec){
    for(Particle& part: particles){
        if(rec.inObject(part.pos)){
            if(world.inBounds(part.pos + dx)){
                part.pos+=dx;
                // std::cout << "moved:!\n";
            }      
        }
    }
};
void Species::move(type_calc3 dx){
    dx = dx.elWiseMult(world.getL());
    for(Particle& part: particles){
        part.pos+=dx;
        if(!world.inBounds(part.pos)){
            part.pos-= dx;
            // std::cout << "moved:!\n";
        }      
    }
    // int np = particles.size();
    // for(int p = 0; p < np; p ++){
    //     if(particles[p].macro_weight == 0){
    //         particles[p] = std::move(particles[np-1]);
    //         p--;
    //         np--;
    //     }
    // }
    // particles.erase(particles.begin() + np, particles.end());

};

void Species::move(type_calc lambda, type_calc A){
    type_calc k = 2*Const::pi/lambda;
    for(Particle& part: particles){
        type_calc dx = A*std::sin(k*part.pos[0])*world.getL()[0];
        part.pos[0]+=dx;
        if(!world.inBounds(part.pos)){
            part.pos[0] -= dx;
            // std::cout << "moved:!\n";
        }      
    }
};


/*for diagnostics and outputs*/
type_calc Species::getMicroCount(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc micro_count{};
    for(Particle& part : particles){
        micro_count += part.macro_weight;
    }
    return micro_count;
};
type_calc3 Species::getMomentum(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc3 mom{};
    for(Particle& part : particles){
        mom += part.macro_weight * part.vel;
    }
    return mom * mass;
};
type_calc Species::getKE(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc ke{};
    for(Particle& part : particles){
        ke += part.macro_weight * (part.vel * part.vel); //using vector multiplication
    }
    return ke * mass * 0.5;
};
type_calc Species::getPE(){ //inmplement somewhere else variable that tracks it so we dont count it every time?
    type_calc pe{};
    type_calc phi_part{};
    type_calc3 lc{};
    for(Particle& part : particles){
        lc = world.XtoL(part.pos);
        phi_part = world.phi.gather(lc);
        pe += part.macro_weight * phi_part;
    }
    return pe * charge;
};
void Species::updateAverages(){
    den_avg.updateMovingAverage(den);
};
void Species::sampleMoments(){
    for(Particle& part: particles){ 
        type_calc3 lc = world.XtoL(part.pos);
        n_sum.scatter(lc, part.macro_weight); //WHY not clearing before?
        nv_sum.scatter(lc, part.macro_weight*part.vel); //TODO check if type_calc -> data_type in scatter works correctly?????
        nuu_sum.scatter(lc, part.macro_weight*part.vel[0]*part.vel[0]);
        nvv_sum.scatter(lc, part.macro_weight*part.vel[1]*part.vel[1]);
        nww_sum.scatter(lc, part.macro_weight*part.vel[2]*part.vel[2]);
    }
};
void Species::computeGasProperties(){
    vel = nv_sum/n_sum;

    for(int i = 0; i < world.ni; i++){
        for(int j = 0; j < world.nj; j++){
            for(int k = 0; k < world.nk; k++){
                type_calc count = n_sum[i][j][k];
                if(count<=0){
                    T[i][j][k] = 0;
                    continue;
                }
                type_calc u_ave = vel[i][j][k][0];
                type_calc v_ave = vel[i][j][k][1];
                type_calc w_ave = vel[i][j][k][2];
                type_calc u2_ave = nuu_sum[i][j][k]/count;
                type_calc v2_ave = nvv_sum[i][j][k]/count;
                type_calc w2_ave = nww_sum[i][j][k]/count;

                type_calc uu = u2_ave - u_ave*u_ave;
                type_calc vv = v2_ave - v_ave*v_ave;
                type_calc ww = w2_ave - w_ave*w_ave;

                T[i][j][k] = mass/(2*Const::k)*(uu + vv + ww);

            }
        }
    }
};
void Species::clearSamples(){
    n_sum.clear();
    nv_sum.clear();
    nuu_sum.clear();
    nvv_sum.clear();
    nww_sum.clear();

};
void Species::computeMacroParticlesCount(){
    macro_part_count.clear();
    for(Particle& part: particles){
        int3 ijk = world.XtoIJK(part.pos);
        macro_part_count[ijk[0]][ijk[1]][ijk[2]] += 1;
    }
};
const std::vector<Particle>& Species::getConstPartRef() const{
    return particles;
};
std::vector<Particle>& Species::getPartRef(){
    return particles;
}; 
const Particle& Species::getConstPartRef(int i) const{
    return particles[i];
};





/*for advance*/
type_calc3 Species::sampleReflectedVelocity(const type_calc3& pos, const type_calc v_mag1, const type_calc3& n) const{
    type_calc v_th = sampleVth(300); //pass Temperature, assuming here 300
    const type_calc a_th = 1; //thermal accommodation coefficient
    type_calc v_mag2 = v_mag1 + a_th*(v_th - v_mag1);

    // part which in book is sphereDiffuseVector
    type_calc sin_theta = rnd();
    type_calc cos_theta = std::sqrt(1 - sin_theta*sin_theta);
    type_calc psi = 2*Const::pi*rnd();

    type_calc3 t1;
    if(n*type_calc3{1,0,0} != 0) t1 = n.cross({1,0,0});
    else t1 = n.cross({0,1,0});
    type_calc3 t2 = n.cross(t1);
    type_calc3 diffuse_vector = sin_theta*std::cos(psi)*t1 + std::sin(psi)*t2 + cos_theta*n;
    // end of that part

    return v_mag2*diffuse_vector;
};

type_calc Species::sampleVth(const type_calc T) const{
    type_calc v_th = sqrt(2*Const::k*T/mass);
    type_calc v1 = v_th * (rnd() + rnd() + rnd() - 1.5);
    type_calc v2 = v_th * (rnd() + rnd() + rnd() - 1.5);
    type_calc v3 = v_th * (rnd() + rnd() + rnd() - 1.5);
    return std::sqrt(v1*v1 + v2*v2 + v3*v3);
};
type_calc3 Species::sampleV3th(const type_calc T) const{
    type_calc v_th = sampleVth(T);
    type_calc theta = 2*Const::pi*rnd(); // random isotropic direction
    type_calc r = -1.0 + 2*rnd();
    type_calc a = std::sqrt(1 - r*r);

    return {v_th*r, v_th*(cos(theta)*a), v_th*(sin(theta)*a)};
};

void Species::sampleVthVariableMpw(const type_calc T, type_calc& set_vel, type_calc& set_mpw) const{
    type_calc v_th = sqrt(2*Const::k*T/mass);

    type_calc v_min = -6 * v_th;
    static const int NUM_BINS = 41;
    static const int PARTS_PER_BIN = 10;
    const int np_uniform = 200;
    const type_calc dv = -2*v_min/NUM_BINS;
    const type_calc a = 1/(std::sqrt(Const::pi) * v_th);

    int bin = (int)(rnd() * NUM_BINS); // values: 0 - 40
    type_calc v = v_min + bin*dv;
    type_calc fm = a*std::exp(-v*v / (v_th*v_th));
    type_calc mpw = (int)(this->mpw0 * np_uniform * fm * dv / PARTS_PER_BIN);

    std::cout << "mpw: " << mpw << " mpw0: " << mpw0 << " bin: " << bin;

    v = v_min + bin*dv + rnd()*dv;

    std::cout << " v: " << v << " v_th: " << v_th << std::endl;

    set_vel = v;
    set_mpw = mpw;  
};

//for sorting
std::vector<std::vector<Particle*>> Species::sort_pointers(){
    std::vector<std::vector<Particle*>> parts_in_cell(world.num_cells);
    for(Particle& part: particles){
        int c = world.XtoC(part.pos);
        parts_in_cell[c].push_back(&part);
    }
    return parts_in_cell;
};
std::vector<std::vector<int>> Species::sort_indexes(){
    std::vector<std::vector<int>> parts_in_cell(world.num_cells);
    for(int p = 0; p < world.num_cells; p++){
        int c = world.XtoC(particles[p].pos);
        parts_in_cell[c].push_back(p);
    }
    return parts_in_cell;
};
void Species::sort_indexes(std::vector<std::vector<int>>& parts_in_cell){

    for(int p = 0; p < world.num_cells; p++){
        int c = world.XtoC(particles[p].pos);
        parts_in_cell[c].push_back(p);
    }
};

// std::unordered_map<int, std::vector<int>> Species::map_indexes(){
//     std::unordered_map<int, std::vector<int>> map_particles_in_cell(world.num_cells);

//     for(int c=0; c < world.num_cells; c++){ //is it worth for performance?
//         int3 ijk = world.CtoIJK(c);
//         int num_parts_in_cell = macro_part_count[ijk[0]][ijk[1]][ijk[2]];
//         if(num_parts_in_cell > 0){
//             map_particles_in_cell[c].reserve(num_parts_in_cell);
//         }
//     }

//     for(int p = 0; p < particles.size(); p++){
//         int c = world.XtoC(particles[c].pos);
//         map_particles_in_cell[c].push_back(p);
//     }
//     return map_particles_in_cell;
// };
void Species::map_indexes(std::unordered_map<int, std::vector<int>>& map_particles_in_cell){ //efficient if most cells are empty

    for(int c=0; c < world.num_cells; c++){ //is it worth for performance? it is with next sorting so there are only rare maps
        int3 ijk = world.CtoIJK(c);
        int num_parts_in_cell = macro_part_count[ijk[0]][ijk[1]][ijk[2]];
        if(!num_parts_in_cell) continue;
        if(num_parts_in_cell > map_particles_in_cell[c].capacity()){
            map_particles_in_cell[c].reserve(num_parts_in_cell); //reserve estimated number of particles (equal if no source works)
        }
    }

    for(int p = 0; p < particles.size(); p++){
        int c = world.XtoC(particles[p].pos);
        map_particles_in_cell[c].push_back(p);
    }
};
void Species::map_indexes(std::unordered_map<int, std::vector<int>>& map_new_particles_in_cell, std::unordered_map<int, std::vector<int>>& already_mapped_other_species_particles){
    for(const std::pair<int, std::vector<int>> pair: already_mapped_other_species_particles){
        int c = pair.first;
        int3 ijk = world.CtoIJK(c);
        int num_parts_in_cell = macro_part_count[ijk[0]][ijk[1]][ijk[2]];
        if(!num_parts_in_cell){
            //already_mapped_other_species_particles.erase(c); cannot do it during 
            continue;
        }
        if(num_parts_in_cell > map_new_particles_in_cell[c].capacity()){
            map_new_particles_in_cell[c].reserve(num_parts_in_cell); //reserve estimated number of particles (equal if no source works)
        }
    }
    auto end = map_new_particles_in_cell.end();
    for(int p = 0; p < particles.size(); p++){
        int c = world.XtoC(particles[p].pos);
        if(map_new_particles_in_cell.find(c)!=end){
            map_new_particles_in_cell[c].push_back(p);
        }
    }
};

