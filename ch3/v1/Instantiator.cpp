#include "Instantiator.h"



Instantiator::Instantiator(std::unique_ptr<World>& otherWorld, std::vector<Species>& otherSpecies, std::unique_ptr<PotentialSolver>& otherSolver){
    std::cout << "WORLD\n";
    readParametersWorld();

    world = std::make_unique<World>(Wnn);
    world->setExtents(Wx0, Wxm);
    world->setTime(Wdt, Wts);
    otherWorld = std::make_unique<World>(std::move(*world));
    std::cout << "SPECIES\n";
    readParametersSpecies();
    otherSpecies.reserve(Pids.size());
    for(int i = 0; i < Pids.size(); i++){
        otherSpecies.emplace_back(Pids[i], Pmasses[i], Pcharges[i], *otherWorld, TEMP_MPW0);
        if(QS){
            otherSpecies[i].loadParticleBoxQS(Pstarts[i], Pends[i], Pdensities[i], Pgrids[i]);
        }
        else{
            otherSpecies[i].loadParticleBox(Pstarts[i], Pends[i], Pdensities[i], Pgrids[i][0]);
        }
    }

    std::cout << "SOLVER\n";
    readParametersSolver();
    otherSolver = std::make_unique<PotentialSolver>(*otherWorld, Smax_it, Stolerance, Stype);
};

void Instantiator::readParametersWorld(){
    std::string task = "read";
    bool first_iteration = true;
    while(task != "n" && task != "next"){
        //pytanie o typ zadania
        if(first_iteration){
            std::cout << "w -> write parameters, d -> use deafult parameters:\n";
        }
        else{
            std::cout << "w -> write all parameters again, c -> correct one parameter, d -> use deafult parameters, n -> next, end writing parameters: ";
        }
        read_value_and_check_type(task);
        lower(task);  

        //opcje dla konkretnych zadan
        if(!task.compare("w") || !task.compare("write")){

            std::cout << Wnames[0] << ": ";
            read_value_and_check_type(Wnn[0], 0, false);
            read_value_and_check_type(Wnn[1], 0, false);
            read_value_and_check_type(Wnn[2], 0, false);
            std::cout << Wnames[1] << ": ";
            read_value_and_check_type(Wx0[0]);
            read_value_and_check_type(Wx0[1]);
            read_value_and_check_type(Wx0[2]);
            std::cout << Wnames[2] << ": ";
            read_value_and_check_type(Wxm[0]);
            read_value_and_check_type(Wxm[1]);
            read_value_and_check_type(Wxm[2]);
            std::cout << Wnames[3] << ": ";
            read_value_and_check_type(Wdt, 0, false);
            std::cout << Wnames[4] << ": ";
            read_value_and_check_type(Wts, 0, false);

            std::cout << Wnames[0] << " = " << Wnn << " "  << Wnames[0] << " = " << Wx0 << " "  << Wnames[0] << " = " << Wxm << " "  << Wnames[0] << " = " << Wdt << " " << Wnames[0] << " = " << Wts <<"\n";


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("correct") || !task.compare("c")){
            for(std::string& n : Wnames){
                std::cout << n << " ";
            }
            std::cout<<"\n";
            std::cout << "Which value to correct? write its name: ";
            std::string name;
            read_value_and_check_type(name);
            lower(name);
            
            if(!name.compare(Wnames[0])){
                std::cout << "wrtie new value for ";
                std::cout << Wnames[0] << ": ";
                read_value_and_check_type(Wnn[0], 0, false);
                read_value_and_check_type(Wnn[1], 0, false);
                read_value_and_check_type(Wnn[2], 0, false);
            }else if(!name.compare(Wnames[1])){
                std::cout << "wrtie new value for ";
                std::cout << Wnames[1] << ": ";
                read_value_and_check_type(Wx0[0]);
                read_value_and_check_type(Wx0[1]);
                read_value_and_check_type(Wx0[2]);
            }else if(!name.compare(Wnames[2])){
                std::cout << "wrtie new value for ";
                std::cout << Wnames[2] << ": ";
                read_value_and_check_type(Wxm[0]);
                read_value_and_check_type(Wxm[1]);
                read_value_and_check_type(Wxm[2]);
            }else if(!name.compare(Wnames[3])){
                std::cout << "wrtie new value for ";
                std::cout << Wnames[3] << ": ";
                read_value_and_check_type(Wdt, 0, false);
            }else if(!name.compare(Wnames[4])){
                std::cout << "wrtie new value for ";
                std::cout << Wnames[4] << ": ";
                read_value_and_check_type(Wts, 0, false);
            }
            
            std::cout << Wnames[0] << " = " << Wnn << " "  << Wnames[0] << " = " << Wx0 << " "  << Wnames[0] << " = " << Wxm << " "  << Wnames[0] << " = " << Wdt << " " << Wnames[0] << " = " << Wts <<"\n";

            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("d") || !task.compare("default")){
            Wnn = {21, 21, 21};
            Wx0 = {-0.1, -0.1, 0.0};
            Wxm = {0.1, 0.1, 0.2};
            Wdt = 2e-10;
            Wts = 1e4;         
            std::cout << Wnames[0] << " = " << Wnn << " "  << Wnames[1] << " = " << Wx0 << " "  << Wnames[2] << " = " << Wxm << " "  << Wnames[3] << " = " << Wdt << " " << Wnames[4] << " = " << Wts <<"\n";


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        if(task.compare("n") && task.compare("next")){
            std::cout << "invalid option, tey again\n";
        }
    }
}

void Instantiator::readParametersSolver(){
    std::string task = "read";
    bool first_iteration = true;
    while(task != "n" && task != "next"){
        //pytanie o typ zadania
        if(first_iteration){
            std::cout << "w -> write parameters, d -> use deafult parameters:\n";
        }
        else{
            std::cout << "w -> write all parameters again, c -> correct one parameter, d -> use deafult parameters, n -> next, end writing parameters: ";
        }
        read_value_and_check_type(task);
        lower(task);  

        //opcje dla konkretnych zadan
        if(!task.compare("w") || !task.compare("write")){

            std::cout << Snames[0] << ": ";
            read_value_and_check_type(Smax_it, 0, false);
            std::cout << Snames[1] << ": ";
            read_value_and_check_type(Stolerance, 0, false);

            std::cout << Snames[0] << " = " << Smax_it << " "  << Snames[1] << " = " << Stolerance << "\n";


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("correct") || !task.compare("c")){
            for(std::string& n : Snames){
                std::cout << n << " ";
            }
            std::cout<<"\n";
            std::cout << "Which value to correct? write its name: ";
            std::string name;
            read_value_and_check_type(name);
            lower(name);
            
            if(!name.compare(Snames[0])){
                std::cout << "wrtie new value for ";
                std::cout << Snames[0] << ": ";            
                read_value_and_check_type(Smax_it, 0, false);
            }else if(!name.compare(Snames[1])){
                std::cout << "wrtie new value for ";
                std::cout << Snames[1] << ": ";
                read_value_and_check_type(Stolerance, 0, false);
            }
            
            std::cout << Snames[0] << " = " << Smax_it << " "  << Snames[1] << " = " << Stolerance <<"\n";

            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("d") || !task.compare("default")){
            Smax_it = 1e4;
            Stolerance = 1e-4;       
            std::cout << Snames[0] << " = " << Smax_it << " "  << Snames[1] << " = " << Stolerance <<"\n";


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        if(task.compare("n") && task.compare("next")){
            std::cout << "invalid option, tey again\n";
        }
    }
}

void Instantiator::readParametersSpecies(){
    std::string task = "read";
    bool first_iteration = true;
    while(task != "n" && task != "next"){
        //pytanie o typ zadania
        if(first_iteration){
            std::cout << "w -> write parameters, d -> use deafult parameters:\n";
        }
        else{
            std::cout << "w -> write all parameters again, c -> correct one parameter, d -> use deafult parameters, n -> next, end writing parameters: ";
        }
        read_value_and_check_type(task);
        lower(task);  

        //opcje dla konkretnych zadan
        if(!task.compare("w") || !task.compare("write")){

            std::string Pid;
            type_calc Pmass;
            type_calc Pcharge;
            int3 Pgrid;
            type_calc3 Pstart;
            type_calc3 Pend;
            type_calc Pdensity;  
            // mass, charge, name
            std::cout << Pnames[0] << ": ";
            read_value_and_check_type(Pid);
            std::cout << Pnames[1] << ": ";
            read_value_and_check_type(Pmass, 0, false);
            std::cout << Pnames[2] << ": ";
            read_value_and_check_type(Pcharge, 0, false);
            std::cout << Pnames[3] << ": ";
            read_value_and_check_type(Pgrid[0], 0, false);
            read_value_and_check_type(Pgrid[1], 0, false);
            read_value_and_check_type(Pgrid[2], 0, false);
            std::cout << Pnames[4] << ": ";
            read_value_and_check_type(Pstart[0]);
            read_value_and_check_type(Pstart[1]);
            read_value_and_check_type(Pstart[2]);
            std::cout << Pnames[5] << ": ";
            read_value_and_check_type(Pend[0]);
            read_value_and_check_type(Pend[1]);
            read_value_and_check_type(Pend[2]);
            std::cout << Pnames[6] << ": ";
            read_value_and_check_type(Pdensity, 0, false);

            std::cout << Pnames[0] << " = " << Pid << " "  << Pnames[1] << " = " << Pmass << " "  << Pnames[2] << " = " << Pcharge << " "  << Pnames[3] << " = " << Pgrid << " " << Pnames[4] << " = " << Pstart << " " << Pnames[5] << " = " << Pend << " " << Pnames[6] << " = " << Pdensity <<"\n";

            Pids.push_back(Pid);
            Pmasses.push_back(Pmass);
            Pcharges.push_back(Pcharge);
            Pgrids.push_back(Pgrid);
            Pstarts.push_back(Pstart);
            Pends.push_back(Pend);
            Pdensities.push_back(Pdensity);
            

            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("correct") || !task.compare("c")){
            if(Pids.empty()){
                std::cout << "nothing to correct\n";
                continue;
            }
            int index_id = 0;
            for(std::string& id: Pids){
                std::cout << index_id << ": " << id << " ";
                index_id++;
            }
            std::cout << "\nwhich element to correct?: ";
            read_value_and_check_type(index_id, 0, (int)Pids.size());
            
            for(std::string& n : Pnames){
                std::cout << n << " ";
            }
            std::cout << "\nwhich value to correct? write it's name: ";

            std::string name;
            read_value_and_check_type(name);
            lower(name);
            if(!name.compare(Pnames[0])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[0] << ": ";            
                read_value_and_check_type(Pids[index_id]);
            }else if(!name.compare(Pnames[1])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[1] << ": ";            
                read_value_and_check_type(Pmasses[index_id], 0, false);
            }else if(!name.compare(Pnames[2])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[2] << ": ";            
                read_value_and_check_type(Pcharges[index_id], 0, false);
            }else if(!name.compare(Pnames[3])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[3] << ": ";            
                read_value_and_check_type(Pgrids[index_id][0], 0, false);
                read_value_and_check_type(Pgrids[index_id][1], 0, false);
                read_value_and_check_type(Pgrids[index_id][2], 0, false);
            }else if(!name.compare(Pnames[4])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[4] << ": ";            
                read_value_and_check_type(Pstarts[index_id][0]);
                read_value_and_check_type(Pstarts[index_id][1]);
                read_value_and_check_type(Pstarts[index_id][2]);
            }else if(!name.compare(Pnames[5])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[5] << ": ";            
                read_value_and_check_type(Pends[index_id][0]);
                read_value_and_check_type(Pends[index_id][1]);
                read_value_and_check_type(Pends[index_id][2]);
            }else if(!name.compare(Pnames[6])){
                std::cout << "wrtie new value for ";
                std::cout << Pnames[6] << ": ";            
                read_value_and_check_type(Pdensities[index_id], 0, false);
            }
            
            std::cout << Pnames[0] << " = " << Pids[index_id] << " "  << Pnames[1] << " = " << Pmasses[index_id] << " "  << Pnames[2] 
            << " = " << Pcharges[index_id] << " "  << Pnames[3] << " = " << Pgrids[index_id] << " " << Pnames[4] << " = " << Pstarts[index_id] 
            << " " << Pnames[5] << " = " << Pends[index_id] << " " << Pnames[6] << " = " << Pdensities[index_id] <<"\n";


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("d") || !task.compare("default")){
            std::string name;
            std::cout << "which element to add? 1= O+, 2= e-: ";
            read_value_and_check_type(name);
            if(!name.compare("1")){
                Pids.push_back("O+");
                Pmasses.push_back(16*Const::amu);
                Pcharges.push_back(Const::q_e);
                Pgrids.emplace_back(41, 41, 41);
                Pstarts.emplace_back(world->getX0());
                Pends.emplace_back(world->getXm());
                Pdensities.push_back(1e11);
            }else if(!name.compare("2")){
                Pids.push_back("e-");
                Pmasses.push_back(Const::m_e);
                Pcharges.push_back(-Const::q_e);
                Pgrids.emplace_back(21, 21, 21);
                Pstarts.emplace_back(world->getX0());
                Pends.emplace_back(world->getXc());
                Pdensities.push_back(1e11); 
            }             


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }else if(!task.compare("d") || !task.compare("del")){
            if(Pids.empty()){
                std::cout << "nothing to delete" << std::endl;
                continue;
            }
            std::cout << "deleted " << Pids.back() << "\n";
            Pids.pop_back();
            Pmasses.pop_back();
            Pcharges.pop_back();
            Pgrids.pop_back();
            Pstarts.pop_back();
            Pends.pop_back();
            Pdensities.pop_back();


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        if(task.compare("n") && task.compare("next")){
            std::cout << "invalid option, tey again\n";
        }
    }
}