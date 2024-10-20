#ifndef FUNKC_H
#define FUNKC_H
#include <vector>
#include <string>
#include <cctype>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <exception>
#include <limits>
#include <memory>
#include <map>
#include <regex>


#include "World.h"
#include "PotentialSolver.h"
#include "Species.h"

std::map<std::string, double> constants = {
    {"amu", Const::amu},
    {"q_e", Const::q_e},
    {"m_e", Const::m_e},
    {"eps_0", Const::eps_0},
    {"k", Const::k},
    {"pi", Const::pi},
    {"eV2K", Const::eV_to_K}
    // Add more constants if needed
};

double parseExpression(const std::string& expression) {
    std::regex regex(R"((\d+(?:\.\d+)?)[\s]*\*[\s]*(\w+))"); // To recognise expressioins like "16*amu"
    std::smatch match;
    size_t idx = 0;
    if (std::regex_match(expression, match, regex)) {
        double number = std::stod(match[1].str(), &idx);
        if(idx != match[1].str().length()){
            throw std::invalid_argument("Podano niewlasciwy typ danych. ");
        }
        std::string constantName = match[2].str();

        // Sprawdź, czy stała istnieje w mapie
        if (constants.find(constantName) != constants.end()) {
            return number * constants[constantName];
        } else {
            throw std::invalid_argument("Nieznana stała: " + constantName);
        }
    }
    else{
        double number = std::stod(expression, &idx);
        if(idx != expression.length()){
            throw std::invalid_argument("Podano niewlasciwy typ danych. ");
        }
        return number;
    }
        //throw std::invalid_argument("Nieprawidłowy format wyrażenia: " + expression);
}

template <class T>
int oom(T x){ // to get order of magnitude
    int temp = x;
    int oom = 0;
    while(temp > 0){
        temp/=10;
        oom++;
    }
    return oom;
};

std::string lower(std::string& str){ //to lower every character in string
    std::transform(str.begin(), str.end(), str.begin(),
       [](unsigned char c){ return std::tolower(c); });
    return str;
};

template <class T, class V, class U>
void read_value_and_check_type(T& value, V min = std::numeric_limits<T>::min(), U max = std::numeric_limits<T>::max(), bool can_be_equal = true){
    max = T(max);
    min = T(min);
    bool wrong_type_flag = true;
    while(wrong_type_flag){
        wrong_type_flag = false;
        std::string input;
        try{
            std::cin >> input;
            value = static_cast<T>(parseExpression(input));

            if(std::cin.fail()) {
                throw std::invalid_argument("Podano niewlasciwy typ danych. ");
            }
            if(value < min) {
                std::stringstream ss;
                ss << "Podano zbyt mala wartosc. min = " << min << " ";
                throw std::out_of_range(ss.str());
            }
            else if(value > max){
                std::stringstream ss;
                ss << "Podano zbyt duza wartosc. max = " << max << " ";
                throw std::out_of_range(ss.str());
            }
            else if(!can_be_equal){
                if(value == min) {
                    std::stringstream ss;
                    ss << "Podano dolna niedozwolona wartosc graniczna. min = " << min << " ";
                    throw std::out_of_range(ss.str());
                }
                else if(value == max){
                    std::stringstream ss;
                    ss << "Podano gorna niedozwolona wartosc graniczna. max = " << max << " ";
                    throw std::out_of_range(ss.str());
                } 
            }
        }catch(const std::invalid_argument& e){
            std::cout << e.what() << "Sproboj ponownie: ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignoruje błędne dane
            wrong_type_flag = true;
        }catch (const std::out_of_range& e) {
            std::cout << e.what() << "Sprobuj ponownie: ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            wrong_type_flag = true;
        }
    }
    std::cin.clear();
}

template <class T, class V>
void read_value_and_check_type(T& value, V min = std::numeric_limits<T>::min(), bool can_be_equal = true){
    read_value_and_check_type(value, min, std::numeric_limits<T>::max(), can_be_equal);
}

template <class T> //used for strings or not cheching boundries
void read_value_and_check_type(T& value){
    bool wrong_type_flag = true;
    while(wrong_type_flag){
        wrong_type_flag = false;
        std::string input;
        try{
            std::cin >> input;
            value = static_cast<T>(parseExpression(input));
            
            if(std::cin.fail()) {
                throw std::invalid_argument("Podano niewlasciwy typ danych. ");
            }

        }catch(const std::invalid_argument& e){
            std::cout << e.what() << "Sproboj ponownie: ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            wrong_type_flag = true;
        }
    }
    std::cin.clear();
}

void read_value_and_check_type(std::string& value){
    bool wrong_type_flag = true;
    while(wrong_type_flag){
        wrong_type_flag = false;
        try{
            std::cin >> value;

            if (std::cin.fail()) {
                throw std::invalid_argument("Podano niewlasciwy typ danych.");
            }

        }catch(const std::invalid_argument& e){
            std::cout << e.what() << "Sproboj ponownie: ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            wrong_type_flag = true;
        }
    }
    std::cin.clear();
}




class Instantiator{
private:
    //world 
    int3                     Wnn    = {21, 21, 21};
    type_calc3               Wx0    = {-0.1, -0.1, 0.0};
    type_calc3               Wxm    = {0.1, 0.1, 0.2};
    type_calc                Wdt    = 2e-10;
    int                      Wts    = 1e4;
    std::vector<std::string> Wnames = {"grid", "x0", "xm", "dt", "ts"};

    //species and particles
    std::vector<std::string> Pids;
    std::vector<type_calc> Pmasses;
    std::vector<type_calc> Pcharges;
    std::vector<int3> Pgrids;
    std::vector<type_calc3> Pstarts;
    std::vector<type_calc3> Pends;
    std::vector<type_calc> Pdensities;
    bool QS = true;
    std::vector<std::string> Pnames = {"n", "m", "cha", "g", "x0", "xm", "d"};

    //solver
    std::vector<std::string> Snames     = {"mi", "tol"};


public:
    //solver 
    /*public bcoz Solver needs a reference to world constructed earlier, so we cant use anything else than constructor with
      reference to world object0, i might have confused reference with something else but point still stands */
    int       Smax_it    = 1e4;
    type_calc Stolerance = 1e-4;
    
    std::unique_ptr<World> world; //needen in readFunctions as a values storage xD 
    //std::vector<Species> species; /not needed
    //std::unique_ptr<PotentialSolver> solver; //not needed

    Instantiator(std::unique_ptr<World>& otherWorld, std::vector<Species>& otherSpecies, std::unique_ptr<PotentialSolver>& otherSolver);
    void readParametersWorld();
    void readParametersSpecies();
    void readParametersSolver();
};

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
        otherSpecies.emplace_back(Pids[i], Pmasses[i], Pcharges[i], *otherWorld);
        if(QS){
            otherSpecies[i].loadParticleBoxQS(Pstarts[i], Pends[i], Pdensities[i], Pgrids[i]);
        }
        else{
            otherSpecies[i].loadParticleBox(Pstarts[i], Pends[i], Pdensities[i], Pgrids[i][0]);
        }
    }

    std::cout << "SOLVER\n";
    readParametersSolver();
    otherSolver = std::make_unique<PotentialSolver>(*otherWorld, Smax_it, Stolerance);
};

void Instantiator::readParametersWorld(){
    std::string task = "read";
    bool first_iteration = true;
    while(task != "n" && task != "next"){
        //pytanie o typ zadania
        if(first_iteration){
            std::cout << "w -> wprowadz parametry, g -> skorzystaj z przykladowych parametrow:\n";
        }
        else{
            std::cout << "w -> wprowadz wszystkie parametry ponownie, p -> popraw jeden parametr, g -> skorzystaj z przykladowych parametrow, n -> zakoncz wprowadzanie parametrow: ";
        }
        read_value_and_check_type(task);
        lower(task);  

        //opcje dla konkretnych zadan
        if(!task.compare("wprowadz") || !task.compare("w") || !task.compare("write")){

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
        else if(!task.compare("popraw") || !task.compare("p") || !task.compare("pop")){
            for(std::string& n : Wnames){
                std::cout << n << " ";
            }
            std::cout<<"\n";
            std::cout << "ktora wartosc poprawic? wprowadz jej nazwe: ";
            std::string name;
            read_value_and_check_type(name);
            lower(name);
            
            if(!name.compare(Wnames[0])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Wnames[0] << ": ";
                read_value_and_check_type(Wnn[0], 0, false);
                read_value_and_check_type(Wnn[1], 0, false);
                read_value_and_check_type(Wnn[2], 0, false);
            }else if(!name.compare(Wnames[1])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Wnames[1] << ": ";
                read_value_and_check_type(Wx0[0]);
                read_value_and_check_type(Wx0[1]);
                read_value_and_check_type(Wx0[2]);
            }else if(!name.compare(Wnames[2])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Wnames[2] << ": ";
                read_value_and_check_type(Wxm[0]);
                read_value_and_check_type(Wxm[1]);
                read_value_and_check_type(Wxm[2]);
            }else if(!name.compare(Wnames[3])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Wnames[3] << ": ";
                read_value_and_check_type(Wdt, 0, false);
            }else if(!name.compare(Wnames[4])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Wnames[4] << ": ";
                read_value_and_check_type(Wts, 0, false);
            }
            
            std::cout << Wnames[0] << " = " << Wnn << " "  << Wnames[0] << " = " << Wx0 << " "  << Wnames[0] << " = " << Wxm << " "  << Wnames[0] << " = " << Wdt << " " << Wnames[0] << " = " << Wts <<"\n";

            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("g") || !task.compare("gotowe")){
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
            std::cout << "niepoprawna opcja, sproboj jeszcze raz\n";
        }
    }
}

void Instantiator::readParametersSolver(){
    std::string task = "read";
    bool first_iteration = true;
    while(task != "n" && task != "next"){
        //pytanie o typ zadania
        if(first_iteration){
            std::cout << "w -> wprowadz parametry, g -> skorzystaj z przykladowych parametrow:\n";
        }
        else{
            std::cout << "w -> wprowadz wszystkie parametry ponownie, p -> popraw jeden parametr, g -> skorzystaj z przykladowych parametrow, n -> zakoncz wprowadzanie parametrow: ";
        }
        read_value_and_check_type(task);
        lower(task);  

        //opcje dla konkretnych zadan
        if(!task.compare("wprowadz") || !task.compare("w") || !task.compare("write")){

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
        else if(!task.compare("popraw") || !task.compare("p") || !task.compare("pop")){
            for(std::string& n : Snames){
                std::cout << n << " ";
            }
            std::cout<<"\n";
            std::cout << "ktora wartosc poprawic? wprowadz jej nazwe: ";
            std::string name;
            read_value_and_check_type(name);
            lower(name);
            
            if(!name.compare(Snames[0])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Snames[0] << ": ";            
                read_value_and_check_type(Smax_it, 0, false);
            }else if(!name.compare(Snames[1])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Snames[1] << ": ";
                read_value_and_check_type(Stolerance, 0, false);
            }
            
            std::cout << Snames[0] << " = " << Smax_it << " "  << Snames[1] << " = " << Stolerance <<"\n";

            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        else if(!task.compare("g") || !task.compare("gotowe")){
            Smax_it = 1e4;
            Stolerance = 1e-4;       
            std::cout << Snames[0] << " = " << Smax_it << " "  << Snames[1] << " = " << Stolerance <<"\n";


            if(first_iteration){
                first_iteration = false;
            }
            continue;
        }
        if(task.compare("n") && task.compare("next")){
            std::cout << "niepoprawna opcja, sproboj jeszcze raz\n";
        }
    }
}

void Instantiator::readParametersSpecies(){
    std::string task = "read";
    bool first_iteration = true;
    while(task != "n" && task != "next"){
        //pytanie o typ zadania
        if(first_iteration){
            std::cout << "w -> wprowadz parametry, g -> skorzystaj z przykladowych parametrow:\n";
        }
        else{
            std::cout << "w -> wprowadz wszystkie parametry ponownie, p -> popraw jeden parametr, g -> skorzystaj z przykladowych parametrow, n -> zakoncz wprowadzanie parametrow: ";
        }
        read_value_and_check_type(task);
        lower(task);  

        //opcje dla konkretnych zadan
        if(!task.compare("wprowadz") || !task.compare("w") || !task.compare("write")){

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
        else if(!task.compare("popraw") || !task.compare("p") || !task.compare("pop")){
            if(Pids.empty()){
                std::cout << "nic do poprawy\n";
                continue;
            }
            int index_id = 0;
            for(std::string& id: Pids){
                std::cout << index_id << ": " << id << " ";
                index_id++;
            }
            std::cout << "\nktora czastke poprawic?: ";
            read_value_and_check_type(index_id, 0, (int)Pids.size());
            
            for(std::string& n : Pnames){
                std::cout << n << " ";
            }
            std::cout << "\nktora wartosc poprawic? wprowadz jej nazwe: ";

            std::string name;
            read_value_and_check_type(name);
            lower(name);
            if(!name.compare(Pnames[0])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Pnames[0] << ": ";            
                read_value_and_check_type(Pids[index_id]);
            }else if(!name.compare(Pnames[1])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Pnames[1] << ": ";            
                read_value_and_check_type(Pmasses[index_id], 0, false);
            }else if(!name.compare(Pnames[2])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Pnames[2] << ": ";            
                read_value_and_check_type(Pcharges[index_id], 0, false);
            }else if(!name.compare(Pnames[3])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Pnames[3] << ": ";            
                read_value_and_check_type(Pgrids[index_id][0], 0, false);
                read_value_and_check_type(Pgrids[index_id][1], 0, false);
                read_value_and_check_type(Pgrids[index_id][2], 0, false);
            }else if(!name.compare(Pnames[4])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Pnames[4] << ": ";            
                read_value_and_check_type(Pstarts[index_id][0]);
                read_value_and_check_type(Pstarts[index_id][1]);
                read_value_and_check_type(Pstarts[index_id][2]);
            }else if(!name.compare(Pnames[5])){
                std::cout << "wprowadz nowa wartosc dla ";
                std::cout << Pnames[5] << ": ";            
                read_value_and_check_type(Pends[index_id][0]);
                read_value_and_check_type(Pends[index_id][1]);
                read_value_and_check_type(Pends[index_id][2]);
            }else if(!name.compare(Pnames[6])){
                std::cout << "wprowadz nowa wartosc dla ";
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
        else if(!task.compare("g") || !task.compare("gotowe")){
            std::string name;
            std::cout << "jaki element dodac? 1= O+, 2= e-: ";
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
        }else if(!task.compare("u") || !task.compare("del") || !task.compare("usun")){
            if(Pids.empty()){
                std::cout << "nic do usuniecia" << std::endl;
                continue;
            }
            std::cout << "usunieto " << Pids.back() << "\n";
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
            std::cout << "niepoprawna opcja, sproboj jeszcze raz\n";
        }
    }
}

#endif
