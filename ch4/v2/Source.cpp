#include "Source.h"

Source::Source(Species& species, World& world, type_calc v_drift, type_calc den, std::string inlet_face_name, type_calc area_frac) noexcept:
 sp{species}, world{world}, v_drift{v_drift}, den{den}, area_frac{area_frac}{
    inlet_face_index = inletName2Index(inlet_face_name);
    if(inlet_face_index<0){
        try{
            throw std::invalid_argument("inlet_face_index invlaid, value:" + inlet_face_name);
        }
        catch(const std::invalid_argument& e){
            std::cerr << e.what() << std::endl;
            inlet_face_index = 0;
        }
    }
    // std::cout << "inlet face index: " << inlet_face_index << "\n";
    dx = world.getDx();
    x0 = world.getX0();
    Lx = dx[0] * (world.ni-1);
    Ly = dx[1] * (world.nj-1);
    Lz = dx[2] * (world.nk-1);
    if(inlet_face_index == 0 || inlet_face_index == 1){
        A = Ly*Lz;        
    }else if(inlet_face_index == 2 || inlet_face_index == 3){
        A = Lx*Lz;
    }else{
        A = Lx*Ly;
    }
    num_micro = den*v_drift*A*world.getDt(); //N = n*v*A*dt
    world.addInlet(inlet_face_name);
};

///////////////////////////////// ColdBeamSource ///////////////////////////////////

ColdBeamSource::ColdBeamSource(Species& species, World& world, type_calc v_drift, type_calc den, std::string inlet_face_name, type_calc area_frac) noexcept:
 Source(species, world, v_drift, den, inlet_face_name, area_frac){
    setSampleStrategy();
};

void ColdBeamSource::setSampleStrategy() noexcept{
    switch(inlet_face_index){
        case 0: //x-
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    sp.addParticle({x0[0], x0[1] + rnd()*Ly, x0[2] + rnd()*Lz}, {v_drift, 0, 0}, sp.mpw0);
                    // sp.addParticle({x0[0], x0[1] + (1-area_frac)*Ly/2 + rnd()*Ly*area_frac, x0[2] + (1-area_frac)*Lz/2 + rnd()*Lz*area_frac}, vel, sp.mpw0);

                }
            };
            break;
        case 1: //x+
            Lx -= dx[0];
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());
                for(int i = 0; i < num_macro; i++){
                    sp.addParticle({x0[0]+Lx, x0[1] + rnd()*Ly, x0[2] + rnd()*Lz}, {-v_drift, 0, 0}, sp.mpw0);
                }
            };
            break;
        case 2: //y-
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());
                for(int i = 0; i < num_macro; i++){
                    sp.addParticle({x0[0] + rnd()*Lx, x0[1], x0[2] + rnd()*Lz}, {0, v_drift, 0}, sp.mpw0);
                }
            };
            break;
        case 3: //y+
            Ly -= dx[1];

            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());
                for(int i = 0; i < num_macro; i++){
                    sp.addParticle({x0[0] + rnd()*Lx, x0[1]+Ly, x0[2] + rnd()*Lz}, {0, -v_drift, 0}, sp.mpw0);
                }
            };
            break;
        case 4: //z-
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());
                for(int i = 0; i < num_macro; i++){
                    sp.addParticle({x0[0] + rnd()*Lx, x0[1] + rnd()*Ly, x0[2]}, {0, 0, v_drift}, sp.mpw0);
                }
            };
            break;
        case 5: //z+
            Lz -= dx[2];
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());
                for(int i = 0; i < num_macro; i++){
                    sp.addParticle({x0[0] + rnd()*Lx, x0[1] + rnd()*Ly, x0[2] + Lz}, {0, 0, -v_drift}, sp.mpw0);
                }
            };
            break;
        

    }
};

void ColdBeamSource::sample()const noexcept{
    sample_strategy();
};

///////////////////////////////// ColdBeamSource ///////////////////////////////////

WarmBeamSource::WarmBeamSource(Species& species, World& world, type_calc v_drift, type_calc den, type_calc T, std::string inlet_face_name, type_calc area_frac) noexcept:
 Source(species, world, v_drift, den, inlet_face_name, area_frac), T{T}{
    setSampleStrategy();
};

void WarmBeamSource::setSampleStrategy() noexcept{ //TODO
    switch(inlet_face_index){
        case 0: //x-
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    type_calc3 vel = sp.sampleV3th(T);
                    vel[0] += v_drift;

                    sp.addParticle({x0[0], x0[1] + rnd()*Ly, x0[2] + rnd()*Lz}, vel, sp.mpw0);
                }
            };
            break;
        case 1: //x+
            Lx -= dx[0];

            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    type_calc3 vel = sp.sampleV3th(T);
                    vel[0] -= v_drift;

                    sp.addParticle({x0[0] + Lx, x0[1] + rnd()*Ly, x0[2] + rnd()*Lz}, vel, sp.mpw0);
                }
            };
            break;
        case 2: //y-
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    type_calc3 vel = sp.sampleV3th(T);
                    vel[1] += v_drift;

                    sp.addParticle({x0[0] + rnd()*Lx, x0[1], x0[2] + rnd()*Lz}, vel, sp.mpw0);
                }
            };
            break;
        case 3: //y+
            Ly -= dx[1];

            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    type_calc3 vel = sp.sampleV3th(T);
                    vel[1] -= v_drift;

                    sp.addParticle({x0[0] + rnd()*Lx, x0[1] + Ly, x0[2] + rnd()*Lz}, vel, sp.mpw0);
                }
            };
            break;
        case 4: //z-
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    type_calc3 vel = sp.sampleV3th(T);
                    vel[2] += v_drift;

                    sp.addParticle({x0[0] + rnd()*Lx, x0[1] + rnd()*Ly, x0[2]}, vel, sp.mpw0);
                }
            };
            break;
        case 5: //z+
            Lz -= dx[2];
            sample_strategy = [this](){
                int num_macro = (int)(num_micro/sp.mpw0 + rnd());                
                for(int i = 0; i < num_macro; i++){
                    type_calc3 vel = sp.sampleV3th(T);
                    vel[2] -= v_drift;

                    sp.addParticle({x0[0] + rnd()*Lx, x0[1] + rnd()*Ly, x0[2] + Lz}, vel, sp.mpw0);
                }
            };
            break;
    }
};

void WarmBeamSource::sample()const noexcept{
    sample_strategy();
};

