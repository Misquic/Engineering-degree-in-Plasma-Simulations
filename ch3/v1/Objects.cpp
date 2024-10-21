#include "Objects.h"
#include <map>

/*constructors*/
Objects::Objects(World& world): world{world}{
};
Objects::Objects(const Objects& other) noexcept: world{other.world}{
    objects = other.objects;
};
Objects::Objects(Objects&& other) noexcept: world{other.world}{
    objects = std::move(other.objects);
}; 

/*methods*/
void Objects::advance() const{
    for(const std::shared_ptr<Object>& obj_ptr: objects){
        if(obj_ptr->movable){
            obj_ptr->pos += obj_ptr->vel * world.getDt(); //moves object
            
            //if forces to objects are coded
            //vel += acc / mass * world.getDt();
        }
    }
};

std::ostream& operator<<(std::ostream& out, Objects objects){
    for(const std::shared_ptr<Object>& obj_ptr: objects.objects){
        out << *obj_ptr;
        out << "\n";
    }
    return out;
};

// std::map<std::string, double> constants = {
//     {"amu", Const::amu},
//     {"q_e", Const::q_e},
//     {"m_e", Const::m_e},
//     {"eps_0", Const::eps_0},
//     {"k", Const::k},
//     {"pi", Const::pi},
//     {"eV2K", Const::eV_to_K}
//     // Add more constants if needed
// };


