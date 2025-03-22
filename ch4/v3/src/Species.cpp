#include <algorithm>
#include <thread>
#include "Species.h"
#include "Rnd.h"
#include "Config.h"
#include "ThreadPool.h"


/*Particle constructors*/
Particle::Particle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight) noexcept: pos{x, y, z}, vel{u, v, w}, macro_weight{macro_weight} {
};
Particle::Particle(type_calc3 pos, type_calc3 vel, type_calc macro_weight) noexcept: pos{pos}, vel{vel}, macro_weight{macro_weight} {
};

Particle::Particle(Particle&& other) noexcept: pos(std::move(other.pos)), vel(std::move(other.vel)), macro_weight(other.macro_weight){
};
Particle::Particle(const Particle& other) noexcept: pos(other.pos), vel(other.vel), macro_weight(other.macro_weight){
};

Particle& Particle::operator=(Particle&& other) noexcept{
    if(this!=&other){
        pos = std::move(other.pos);
        vel = std::move(other.vel);
        macro_weight = other.macro_weight;
    }
    return *this;
};

/*Species constructors*/
Species::Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0): name{name}, mass{mass},
 charge{charge}, world{world}, mpw0{mpw0}, den{world.nn}, den_avg{world.nn}, n_sum{world.nn}, nv_sum{world.nn}, nuu_sum{world.nn},
 nvv_sum{world.nn}, nww_sum{world.nn}, T{world.nn}, vel{world.nn}, macro_part_count{world.ni -1, world.nj -1, world.nk-1} {
    sorted_particles_indexes.reserve((world.ni-1)*(world.nj-1)*(world.nk-1));
};
Species::Species(std::string name, type_calc mass, type_calc charge, World& world, type_calc mpw0, type_calc E_ion): name{name}, mass{mass},
 charge{charge}, world{world}, mpw0{mpw0}, E_ion{E_ion}, den{world.nn}, den_avg{world.nn}, n_sum{world.nn}, nv_sum{world.nn}, nuu_sum{world.nn},
 nvv_sum{world.nn}, nww_sum{world.nn}, T{world.nn}, vel{world.nn}, macro_part_count{world.ni -1, world.nj -1, world.nk-1} {
    sorted_particles_indexes.reserve((world.ni-1)*(world.nj-1)*(world.nk-1));
};


/*Species methods*/
size_t Species::getNumParticles() const{
    return particles.size();
};

void Species::advanceNonElectron(Species& neutrals, Species& spherium, type_calc dt){
    Config& config = Config::getInstance();
    // if(config.getMULTITHREADING() && advance_time_multi < advance_time_serial){
    if(false){
        // type_calc time_multi_start = world.getWallTime();
        // if(charge!=0){//ions
        //     if(config.getSPUTTERING()){
        //         advanceSputteringMultithreading(neutrals, spherium, dt);
        //     }else{
        //         advanceNoSputteringMultithreading(neutrals, dt);
        //     }
        // }else{//neutrals
        //     advanceNoSputteringMultithreading(neutrals, dt);
        // }
        // advance_time_multi = world.getWallTime() - time_multi_start;
        // //implement multithreading for sputt and no sput;
    }else{   
        type_calc time_serial_start = world.getWallTime();
        if(charge!=0){//ions
            if(config.getSPUTTERING()){
                advanceSputteringSerial(neutrals, spherium, dt);
            }else{
                advanceNoSputteringSerial(neutrals, dt);
            }
        }else{//neutrals
            advanceNoSputteringSerial(neutrals, dt);
        }
        advance_time_serial = world.getWallTime() - time_serial_start;
    }

};

// void Species::advanceNoSputteringMultithreading()

void Species::advanceSputteringSerial(Species& neutrals, Species& spherium, type_calc dt){//TODO change it so species take into account if its neutral or not? -> subclasses
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
                    type_calc mpw_ratio = part.macro_weight/neutrals.mpw0;

                    /*inject neutrals*/
                    int mp_create = (int)(mpw_ratio + rnd()); //number of macroparticles to create
                    for(int i = 0; i < mp_create; i ++){
                        type_calc3 vel = sampleReflectedVelocity(part.pos, v_mag_pre_impact, n);
                        neutrals.addParticle(part.pos, vel);
                    }
                    if(Config::getInstance().getSPUTTERING()){
                        /*emit sputtered material using simple yield model*/
                        type_calc sput_yield = (v_mag_pre_impact>5e3)?0.1:0;
                        type_calc sput_macroweight_ratio = sput_yield*part.macro_weight/spherium.mpw0;
                        int sput_mp_create = (int)(sput_macroweight_ratio + rnd());//number of sputtered macroparticles to create
                        for(int i = 0; i < sput_mp_create; i++){
                            type_calc3 vel = sampleReflectedVelocity(part.pos, v_mag_pre_impact, n);
                            spherium.addParticle(part.pos, vel);
                        }
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
        #ifdef DEBUG
        if(part.pos.isNan()){
            std::cerr << "Advance Sputtering stop pos nan\n";
        }
        if(part.vel.isNan()){
            std::cerr << "Advance Sputtering stop vel nan\n";
        }
        #endif
    }
    sorted = false;
    particles.erase(particles.begin() + np, particles.end());
};
void Species::advanceNoSputteringSerial(Species& neutrals, type_calc dt){//TODO change it so species take into account if its neutral or not? -> subclasses
    type_calc3 x0 = world.getX0();
    type_calc3 xm = world.getXm();

    type_calc3 lc{};
    type_calc3 ef_part{};
    
    size_t np = particles.size();

    for(size_t p = 0; p < np; p++){
        Particle& part = particles[p];
        #ifdef DEBUG
        if(part.pos.isNan()){
            std::cerr << "start pos nan\n";
        }
        if(part.vel.isNan()){
            std::cerr << "start vel nan\n";
        }
        #endif
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
                    type_calc mpw_ratio = part.macro_weight/neutrals.mpw0;

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
        #ifdef DEBUG
        if(part.pos.isNan()){
            std::cerr << "stop pos nan\n";
        }
        if(part.vel.isNan()){
            std::cerr << "stop vel nan\n";
        }
        #endif
    }
    sorted = false;
    particles.erase(particles.begin() + np, particles.end());
};

void Species::advanceElectrons(type_calc dt){
    type_calc time_start = world.getWallTime();

    Config& config = Config::getInstance();
    size_t num_part_start = particles.size(); 
    if(config.getMULTITHREADING() && advance_time_multi < advance_time_serial){
        dmsg("used multithreading advanceElectron\n");

        type_calc time_multi_start = world.getWallTime();
        // unsigned int num_threads = config.getNumThreads();

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ThreadPool& thread_pool = SingletonThreadPool::getInstance();
        thread_pool.waitForCompletion();
        std::vector<size_t> indexes = splitIntoChunks(particles.size());
        std::vector<std::vector<size_t>> parts_to_delete(indexes.size() - 1);//vector of buffers to delete particles of that indexes
        // std::queue<std::future<void>> results;
        // std::cout << "pre: " << thread_pool.queueSize();
        for(unsigned int i = 0; i < indexes.size()-1; i++){
            // thread_pool.AddTask(operatorHelper, v, ret, indexes[i], indexes[i+1]);
            thread_pool.AddTask([this](unsigned int thread_id, type_calc dt, size_t index_start, size_t index_stop, std::vector<size_t>& t_buff){this->advanceElectronsMultithreading(thread_id, dt, index_start, index_stop, t_buff);},i, dt, indexes[i], indexes[i+1], parts_to_delete[i] );
            // results.emplace(thread_pool.AddTask([this](const tcvector* v_lamb, tcvector* ret_lamb, size_t start, size_t stop){this->operatorHelper(v_lamb, ret_lamb, start, stop);}, &v, &ret, indexes[i], indexes[i+1]));
        }
        // std::cout << "1\n";
        thread_pool.waitForCompletion();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //deleting particles
        size_t number_of_parts_to_delete{};
        for(const std::vector<size_t>& buf: parts_to_delete){
            number_of_parts_to_delete+=buf.size();
        }
        if(number_of_parts_to_delete>0){
            dmsg("number of parts_to delete: " << number_of_parts_to_delete << "\n");
            // std::cout << "number of parts_to delete: " << number_of_parts_to_delete << "\n";
            std::vector<size_t> parts_to_delete_total;
            parts_to_delete_total.reserve(number_of_parts_to_delete);

            for(const std::vector<size_t>& buf: parts_to_delete){//merge buffers of parts to delete
                std::move(buf.begin(), buf.end(), std::back_inserter(parts_to_delete_total));
            }

            std::sort(parts_to_delete_total.begin(), parts_to_delete_total.end(), std::greater<size_t>());//sort descending
            size_t last_index = particles.size()-1;
            for(size_t p: parts_to_delete_total){
                particles[p] = std::move(particles[last_index]);
                last_index--;
            }
            particles.erase(particles.end() - parts_to_delete_total.size(), particles.end());
        }
        advance_time_multi = (world.getWallTime()-time_multi_start)/num_part_start;
    }else{
        dmsg("used serial advanceElectron\n");
        type_calc time_serial_start = world.getWallTime();
        advanceElectronsSerial(dt);
        advance_time_serial = (world.getWallTime()-time_serial_start)/num_part_start;
    }
    dmsg("\ntime electrons [s] : " << world.getWallTime()-time_start);
};
void Species::advanceElectronsMultithreading(unsigned int thread_id, type_calc dt, size_t index_start, size_t index_stop, std::vector<size_t>& parts_to_delete){

    dmsg("advance electrons worker "<< thread_id << "start\n");
    for(size_t p = index_start; p < index_stop; p++){
        Particle& part = particles[p];
#ifdef DEBUG
        for(int i = 0; i < 3; i++){
            if(std::isnan(part.pos[i])){
                std::cerr << "thread_id " << thread_id << " nan\n";
            }
        }
#endif

        type_calc3 lc = world.XtoL(part.pos);
        type_calc3 ef_part = world.ef.gather(lc);

        part.vel += ef_part*(dt*charge/mass);
        part.pos += part.vel*dt;

        int in_object = world.inObject(part.pos);

        if(!world.inBounds(part.pos)){
            parts_to_delete.push_back(p);
        }
        else if(in_object){                    
            parts_to_delete.push_back(p);
        }
#ifdef DEBUG
        for(int i = 0; i < 3; i++){
            if(std::isnan(part.pos[i])){
                std::cerr << "nan\n";
            }
        }
#endif
    }
    dmsg("advance electrons worker "<< thread_id << "finish\n");

};
void Species::advanceElectronsSerial(type_calc dt){
    size_t np = particles.size();

    for(size_t p = 0; p < np; p++){
        Particle& part = particles[p];
        #ifdef DEBUG
        for(int i = 0; i < 3; i++){
            if(std::isnan(part.pos[i])){
                std::cerr << "nan\n";
            }
        }
        #endif

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
        #ifdef DEBUG
        for(int i = 0; i < 3; i++){
            if(std::isnan(part.pos[i])){
                std::cerr << "nan\n";
            }
        }
        #endif
    }
    sorted = false;
    particles.erase(particles.begin() + np, particles.end());
}

void Species::computeNumberDensity(){ //to ommit doind this thice implement using n_sum

    den.clear();
    for(Particle& part : particles)
    {
        type_calc3 lc = world.XtoL(part.pos);
        #ifdef DEBUG
        if(lc.isNan()){
            std::cerr << "pos: " << part.pos << " weigth: " << part.macro_weight << std::endl;
        }
        #endif
        den.scatter(lc, part.macro_weight);
    }

    den/=world.node_vol;
};
void Species::addParticle(type_calc x, type_calc y, type_calc z, type_calc u, type_calc v, type_calc w, type_calc macro_weight){
    addParticle({x, y, z}, {u, v, w}, macro_weight);
};
void Species::addParticle(type_calc3 pos, type_calc3 vel, type_calc macro_weight){
    if(pos.isNan() || vel.isNan()){
        return;
    }
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

    std::cout << "Loading particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weight of macroparticle: " << macro_weight << "\n";

    particles.reserve(num_macro);

    type_calc3 pos;
    type_calc3 vel;
    for(int i = 0; i < num_macro; i++){
        pos = {rnd(x1[0], x2[0]), rnd(x1[1], x2[1]), rnd(x1[2], x2[2])};
        vel = {0, 0, 0}; //stationary particle;

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

    std::cout << "Loading particle " << name << " in a box of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weight of macroparticle: " << this->mpw0 << std::endl;
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

    std::cout << "Loading particle " << name << " in a sphere of volume : "<< box_vol << " with number of microparticles: " << num_micro << " with number of macroparticles: " << num_macro << " and weight of macroparticle: " << this->mpw0 << std::endl;
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
        }      
    }
};
bool Species::setSorted(bool is_sorted) noexcept{
    sorted = is_sorted;
    return sorted;
}
bool Species::isSorted() noexcept{
    return sorted;
}


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

    dmsg("Species::sampleVthVariableMpw:  mpw: " << mpw << " mpw0: " << mpw0 << " bin: " << bin);

    v = v_min + bin*dv + rnd()*dv;

    dmsg(" v: " << v << " v_th: " << v_th << "\n");

    set_vel = v;
    set_mpw = mpw;  
};

//for sorting
std::vector<std::vector<Particle*>> Species::sortPointers(){
    std::vector<std::vector<Particle*>> parts_in_cell(world.num_cells);
    for(Particle& part: particles){
        int c = world.XtoC(part.pos);
        parts_in_cell[c].push_back(&part);
    }
    return parts_in_cell;
};
std::vector<std::vector<int>> Species::sortIndexes(){
    std::cerr << "vector sorting " << name << "\n";
    std::vector<std::vector<int>> parts_in_cell(world.num_cells);
    for(int p = 0; p < particles.size(); p++){
        int c = world.XtoC(particles[p].pos);
        parts_in_cell[c].push_back(p);
    }
    sorted = true;
    sorted_particles_indexes = parts_in_cell;
    return parts_in_cell;
};
std::vector<std::vector<int>> Species::getSortedIndexes(){
    return sorted_particles_indexes;
};

void Species::sortIndexes(std::vector<std::vector<int>>& parts_in_cell){
    std::cerr << "vector sorting " << name << "\n";
    for(int p = 0; p < particles.size(); p++){
        int c = world.XtoC(particles[p].pos);
        parts_in_cell[c].push_back(p);
    }
    sorted = true;
    sorted_particles_indexes = parts_in_cell;

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
void Species::mapIndexes(std::unordered_map<int, std::vector<int>>& map_particles_in_cell){ //efficient if most cells are empty
    std::cerr << "map sorting " << name << "\n";

    for(int p = 0; p < particles.size(); p++){
        int c = world.XtoC(particles[p].pos);
        map_particles_in_cell[c].push_back(p);
    }
};
void Species::mapIndexes(std::unordered_map<int, std::vector<int>>& map_new_particles_in_cell, std::unordered_map<int, std::vector<int>>& already_mapped_other_species_particles){
    std::cerr << "map sorting " << name << "\n";
    
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

std::vector<std::vector<std::vector<std::vector<int>>>> Species::sortVelocitiesInCell(const std::vector<int>& particles_in_cell){
    
    std::vector<std::vector<std::vector<std::vector<int>>>> velocity_grid(m_vel_grid_n[0], std::vector<std::vector<std::vector<int>>>(int(m_vel_grid_n[1]),
     std::vector<std::vector<int>>(int(m_vel_grid_n[2]), std::vector<int>())));
    type_calc3 vel_min = particles[0].vel, vel_max = particles[0].vel;
    
    //find vel min and max
    for(int p: particles_in_cell){
        const type_calc3& vel = particles[p].vel;
        for(int i = 0; i < 3; i++){
            if(vel[i] < vel_min[i]){
                vel_min[i] = vel[i];
            } 
            else if(vel[i] > vel_max[i]){
                vel_max[i] = vel[i];
            }
        }
    }
    // for(int i = 0; i < 3; i++){ // ??
    //     if(vel_max[i] == 0){
    //         std::stringstream ss;
    //         ss << "Species::sort_velocities_in_cell: vel_max[i] == 0, i: " << i << "\n";
    //         throw std::runtime_error(ss.str());
    //     }
    // }
        
    for(int i = 0; i < 3; i++){// to omit situation when velocity of particle = vel_max which results in index = vel_grid_n which overflows array
        if(vel_max[i]>0){
            vel_max[i]*=1.001;
        }else{
            vel_max[i]*=0.999;
        }
    }
    
    type_calc3 vel_range = vel_max-vel_min;
    for(int i = 0; i < 3; i ++){
        if(vel_max[i] == vel_min[i]){
            std::cerr << "Species::sort_velocities_in_cell: vel_min[i] == vel_max[i], i: " << i << "\n"; 
        }
    }
    type_calc3 d_vel = vel_range/m_vel_grid_n;

    for(int p: particles_in_cell){
        type_calc3 ijk = (particles[p].vel - vel_min)/d_vel; 
        for (int i = 0; i < 3; i++)
        {    
            if(ijk[i]>=15)
            {
                std::cout<<"asdfghj\n";
            }
        }
        velocity_grid[int(ijk[0])][int(ijk[1])][int(ijk[2])].push_back(p);
    }
    return velocity_grid;
};

void Species::merge(){
    std::stringstream out;
    out << "Merging " << name << "\nBefore merge: particles size: " << particles.size() << ", energy: " << getKE() << ", momentum: " << getMomentum() << "\n";
    // bool use_map = should_use_map();
    bool use_map = false;
    if(use_map){

    }else{
        if(!sorted){
            sorted_particles_indexes = sortIndexes();
        }

        for(std::vector<int>& cell: sorted_particles_indexes){ // loop through all space cells
            if(cell.size() < 10){
                continue;
            }
            std::vector<std::vector<std::vector<std::vector<int>>>> sorted_velocities_in_cell = sortVelocitiesInCell(cell); //3D vector of vectors to particles indexes

            for(int i = 0; i < m_vel_grid_n[0]; i++){ //loop through all velocity cells
                for(int j = 0; j < m_vel_grid_n[1]; j++){
                    for(int k = 0; k < m_vel_grid_n[2]; k++){
                        if(sorted_velocities_in_cell[i][j][k].size() > 2){ //merge if number of particles in velocity cell > 2
                            std::vector<int>& velocities_cell = sorted_velocities_in_cell[i][j][k];  
                            type_calc total_weight{};
                            type_calc3 total_momentum{}, total_energy{};
                            for(int p: velocities_cell){
                                Particle& particle = particles[p];
                                total_weight += particle.macro_weight;
                                total_momentum += particle.macro_weight*particle.vel;
                                total_energy += particle.macro_weight*(particle.vel.elWiseMult(particle.vel));
                                particle.macro_weight = 0;
                            } 

                            type_calc w_a = total_weight*0.5, w_b = total_weight*0.5; //not always int
                            total_momentum /= total_weight;
                            total_energy /= total_weight;
                            total_energy -= total_momentum.elWiseMult(total_momentum);
                            for(int i = 0; i < 3; i++){
                                if(total_energy[i] < 0){
                                    if(total_energy[i] < -10e-5){
                                        dmsg( "Species::merge: " << name << " total energy < -10e-5 shifted to 0" << total_energy << "\n");
                                        break;
                                    }
                                    total_energy[i] = 0;
                                }
                            }
                            type_calc3 total_energy_sqrt = {std::sqrt(total_energy[0]), std::sqrt(total_energy[1]), std::sqrt(total_energy[2])}; 
                            type_calc3 v_a = total_momentum + total_energy_sqrt;
                            type_calc3 v_b = total_momentum - total_energy_sqrt;

                            // if(total_energy_sqrt.isNan("total_energy_sqrt")){
                            //     type_calc total_weight_temp;
                            //     type_calc3 total_momentum_temp{}, total_energy_temp{};
                            //     for(int p: velocities_cell){
                            //         Particle& particle = particles[p];
                            //         total_weight_temp += particle.macro_weight;
                            //         total_momentum_temp += particle.macro_weight*particle.vel;
                            //         type_calc3 temp2 = particle.macro_weight*(particle.vel.elWiseMult(particle.vel));
                            //         type_calc3 temp = (particle.vel.elWiseMult(particle.vel));
                            //         total_energy_temp += particle.macro_weight*(particle.vel.elWiseMult(particle.vel));
                            //     } 
                            // }
                            // for(int p: velocities_cell){
                            //     Particle& particle = particles[p];
                            //     particle.macro_weight = 0;
                            // } 

                            //sampling position
                            int cell_size = velocities_cell.size();
                            int p_a = velocities_cell[rnd(0, cell_size)];
                            int p_b = velocities_cell[rnd(0, cell_size)];
                            while(p_b == p_a){
                                p_b = velocities_cell[rnd(0, cell_size)];
                            }
                            // std::cerr << "size: " << cell_size << " p_a: " << p_a << " p_b: " << p_b << std::endl;
                            type_calc3 pos_a = particles[p_a].pos;
                            type_calc3 pos_b = particles[p_b].pos;
                            
                            #ifdef DEBUG
                            pos_a.isNan("pos_a");
                            v_a.isNan("v_a");
                            pos_b.isNan("pos_b");
                            v_b.isNan("v_b");
                            #endif
                            addParticle(pos_a, v_a, w_a);
                            addParticle(pos_b, v_b, w_b);
                            //std::cerr << "Merged particle!\n";
                        }
                    }
                }
            }
        }
        int np = particles.size();
        for(int p = 0; p < np; p++){
            if(particles[p].macro_weight == 0){
                particles[p] = std::move(particles[np-1]);
                np--;
                p--;
            }
        }
        if(np != particles.size()){
            sorted = false;
        }
        particles.erase(particles.begin() + np, particles.end());
    }

    out << name << " after merge: particles size: " << particles.size() << ", energy: " << getKE() << ", momentum: " << getMomentum() << "\n";
    std::cout << out.str();
};

#define map_threshold 0.05

bool Species::shouldUseMap(){
    
    int occupied_cells = 0;
    for(int i = 0; i < world.ni_1; i ++){
        for(int j = 0; j < world.nj_1; j ++){
            for(int k = 0; k < world.nk_1; k ++){
                if(macro_part_count[i][j][k]>0){
                    occupied_cells++;
                }
            }
        }
    }
    dmsg(static_cast<type_calc>(occupied_cells)/macro_part_count.size());
    // std::cerr << static_cast<type_calc>(occupied_cells)/macro_part_count.size() << ", ";
    // std::cerr << "\nfraction of occupied cells: " << static_cast<type_calc>(occupied_cells)/macro_part_count.size() << "\n";
    if(static_cast<type_calc>(occupied_cells)/macro_part_count.size() > map_threshold){ //most occupied
        return false;
    }
    return true; //most empty
};


