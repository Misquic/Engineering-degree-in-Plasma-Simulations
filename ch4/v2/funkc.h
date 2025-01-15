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
#include <map>
#include <regex>
#include "all.h"

void print_help();

template <typename T>
T parseArgument(const std::vector<std::string>& args, const std::string& option, T defaultValue) { //to check if command occurs and to read passed value if it is
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == option && i + 1 < args.size()) { //TODO add possibility of more then one option, ends witd something that starts with "--""
            std::istringstream iss(args[i + 1]);
            T value;
            if (iss >> value){
                std::cout << " option: " + option + " = " << value << "\n";
                return value;
            }
        }
    }
    std::cout << " option: " + option + " = " << defaultValue << "\n";
    return defaultValue;
}

bool parseArgument(const std::vector<std::string>& args, const std::string& option);//to chceck if command occurs
double parseExpression(const std::string& expression); // to read constans like eps_0 from string
std::string lower(std::string& str);
void read_value_and_check_type(std::string& value); //for instantiator



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
                throw std::invalid_argument("Passed wrong type. ");
            }
            if(value < min) {
                std::stringstream ss;
                ss << "Passed value too low. min = " << min << " ";
                throw std::out_of_range(ss.str());
            }
            else if(value > max){
                std::stringstream ss;
                ss << "Passed value too high. max = " << max << " ";
                throw std::out_of_range(ss.str());
            }
            else if(!can_be_equal){
                if(value == min) {
                    std::stringstream ss;
                    ss << "Passed value equals forbidden minimal value. min = " << min << " ";
                    throw std::out_of_range(ss.str());
                }
                else if(value == max){
                    std::stringstream ss;
                    ss << "Passed value equals forbidden maximal value. max = " << max << " ";
                    throw std::out_of_range(ss.str());
                } 
            }
        }catch(const std::invalid_argument& e){
            std::cout << e.what() << "Try again: ";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignoruje błędne dane
            wrong_type_flag = true;
        }catch (const std::out_of_range& e) {
            std::cout << e.what() << "Try again: ";
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
                throw std::invalid_argument("Paseed wrong type. ");
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

template<class data_type>
std::ostream& operator<<(std::ostream& out, std::vector<data_type>& vec){
    size_t size = vec.size();
    for(int i = 0; i < size; i ++){
        out << vec[i];
        if(i!=size-1){
            out << ", ";
        }
    }
    return out;
}



#endif
