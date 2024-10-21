#ifndef OBJECTS_H
#define OBJECTS_H
#include <memory>
#include <vector>
#include "Object.h"


//class that is container for Objects

class Objects{
private:
    std::vector<std::shared_ptr<Object>> objects; //container for Object(s) shared pointer so no duplicates of Object(s) are created, we wnat to refer to the same object
    World& world; //reference do world
public:
    Objects(World& world);
    Objects(const Objects& other) noexcept;  //copying constructor
    Objects(Objects&& other) noexcept;       //moving constructor
    ~Objects() noexcept = default;

    /*methods*/
    void advance() const; //for advancing Object(s) in simulation
    template <class T, class... Args>
    void addObject(Args&&... args); //for adding Object(s) to container

    friend std::ostream& operator<<(std::ostream& out, Objects objects);
};

template <class T, class... Args>
void Objects::addObject(Args&&... args){
    static_assert(std::is_base_of<Object, T>::value, "T must derive from Object");
    try{
        objects.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error adding object, supposedly you used wrong arguments that do not match any of constructor of Object or derived class: " << e.what() << '\n';
    } catch (...) {
        std::cerr << "Unknown error occurred while adding object\n";
    }

}




#endif