#include "funkc.h"


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

double parseExpression(const std::string& expression){
    std::regex regex(R"((\d+(?:\.\d+)?)[\s]*\*[\s]*(\w+))"); // To recognise expressioins like "16*amu"
    std::smatch match;
    size_t idx = 0;
    if (std::regex_match(expression, match, regex)) {
        double number = std::stod(match[1].str(), &idx);
        if(idx != match[1].str().length()){
            throw std::invalid_argument("Passed wrong type. ");
        }
        std::string constantName = match[2].str();

        // Sprawdź, czy stała istnieje w mapie
        if (constants.find(constantName) != constants.end()) {
            return number * constants[constantName];
        } else {
            throw std::invalid_argument("Uknown constant: " + constantName);
        }
    }
    else{
        double number = std::stod(expression, &idx);
        if(idx != expression.length()){
            throw std::invalid_argument("Passed wrong type. ");
        }
        return number;
    }
        //throw std::invalid_argument("Nieprawidłowy format wyrażenia: " + expression);
};


std::string lower(std::string& str){ //to lower every character in string
    std::transform(str.begin(), str.end(), str.begin(),
       [](unsigned char c){ return std::tolower(c); });
    return str;
};


void read_value_and_check_type(std::string& value){
    bool wrong_type_flag = true;
    while(wrong_type_flag){
        wrong_type_flag = false;
        try{
            std::cin >> value;

            if (std::cin.fail()) {
                throw std::invalid_argument("Passed wrong type.");
            }

        }catch(const std::invalid_argument& e){
            std::cout << e.what() << "Try again: ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            wrong_type_flag = true;
        }
    }
    std::cin.clear();
}

