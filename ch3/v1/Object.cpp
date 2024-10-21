#include "Object.h"
#include <cmath>

/*constructors*/
Object::Object(type_calc3 pos, type_calc3 vel): pos{pos}, vel{vel}, movable{true}{
};
Object::Object(type_calc3 pos): pos{pos}, vel{vel}, movable{false}{
    vel.clear();
};
Object::Object(const Object& other) noexcept: pos{other.pos}, vel{other.vel}, movable{other.movable}{ //does it work as intended?
};
Object::Object(Object&& other) noexcept: pos{std::move(other.pos)}, vel{std::move(other.vel)}, movable{other.movable}{
};

/*destructors*/
// Object::~Object(){ //not needed here, because type_calc3 is type of Vec3 which has its own destructor
// }

/*operators*/
Object& Object::operator=(const Object& other) noexcept{
    if(this!= &other){
        this->pos = other.pos;
        if(this->movable){
            this->vel = other.vel;
        }
    }
    return *this;
};
Object& Object::operator=(Object&& other) noexcept{
    if(this!= &other){
        this->pos = std::move(other.pos);
        if(this->movable){
            this->vel = std::move(other.vel);
        }
    }
    return *this;
};

/*methods*/
void Object::setPhi(type_calc phi) noexcept{
    this->phi = phi;
};
bool Object::isMovable() noexcept{
    return movable;
}


/*friends*/
std::ostream& operator<<(std::ostream& out, const Object& obj){
    out << "Name: " << std::setw(9) << obj.name;
    out << " pos: " << obj.pos;
    if(obj.movable){
        out << " vel: " << obj.vel;
    }
    return out;
};

//////////////////////////////////////////////////////// SPHERE ////////////////////////////////////////////////////////

/*constructors*/
Sphere::Sphere(type_calc3 pos, type_calc3 vel, type_calc radius): Object(pos, vel), radius{fabs(radius)}{
    name = "Sphere";
};
Sphere::Sphere(type_calc3 pos, type_calc radius): Object(pos), radius{fabs(radius)}{
    name = "Sphere";
};
Sphere::Sphere(const Sphere& other): Object(other), radius{other.radius}{
    name = "Sphere";
};
Sphere::Sphere(Sphere&& other): Object(std::move(other)), radius{other.radius}{
    name = "Sphere";
};

/*operators*/
Sphere& Sphere::operator=(const Object& other) noexcept{
    if(this!= &other){
        Object::operator=(other);
        if(const Sphere* spherePtr = dynamic_cast<const Sphere*>(&other)){
            this->radius = spherePtr->radius;
        }
        // else radius is left as it is
    }
    return *this;
};
Sphere& Sphere::operator=(Object&& other) noexcept{
    if(this!= &other){
        Object::operator=(std::move(other));
        if(const Sphere* spherePtr = dynamic_cast<const Sphere*>(&other)){
            this->radius = spherePtr->radius;
        }
        // else radius is left as it is
    }
    return *this;
};
Sphere& Sphere::operator=(const Sphere& other) noexcept{
    if(this!= &other){
        Object::operator=(other);
        this->radius = other.radius;
    }
    return *this;
};  //uses Object copy operator= and copies radius
Sphere& Sphere::operator=(Sphere&& other) noexcept{
    if(this!= &other){
        Object::operator=(std::move(other));
        this->radius = other.radius;
    }
    return *this;
};       //uses Object move operator= and copies radius

std::ostream& operator<<(std::ostream& out, const Sphere& obj){ 
    out << static_cast<const Object&>(obj);
    out << " radius: " << obj.radius;
    return out;
};

//////////////////////////////////////////////////////// Rectangle ////////////////////////////////////////////////////////

/*constructors*/
Rectangle::Rectangle(type_calc3 pos, type_calc3 vel, type_calc3 sides, type_calc3 orientation): 
 Object(pos, vel), orientation{orientation}{
    this->sides = {fabs(sides[0]), fabs(sides[1]), fabs(sides[2])};
    name = "Rectangle";
};
Rectangle::Rectangle(type_calc3 pos, type_calc3 sides, type_calc3 orientation): 
 Object(pos), orientation{orientation}{
    this->sides = {fabs(sides[0]), fabs(sides[1]), fabs(sides[2])};
    name = "Rectangle";
};
Rectangle::Rectangle(const Rectangle& other): Object(other), sides{other.sides}, orientation{other.orientation}{
    name = "Rectangle";
};
Rectangle::Rectangle(Rectangle&& other): Object(std::move(other)), sides{other.sides}, orientation{other.orientation}{
    name = "Rectangle";
};

/*operators*/
Rectangle& Rectangle::operator=(const Object& other) noexcept{
    if(this!= &other){
        Object::operator=(other);
        if(const Rectangle* RectanglePtr = dynamic_cast<const Rectangle*>(&other)){
            this->sides = RectanglePtr->sides;
            this->orientation = RectanglePtr->orientation;
        }
        // else radius is left as it is
    }
    return *this;
};
Rectangle& Rectangle::operator=(Object&& other) noexcept{
    if(this!= &other){
        Object::operator=(std::move(other));
        if(const Rectangle* RectanglePtr = dynamic_cast<const Rectangle*>(&other)){
            this->sides = RectanglePtr->sides;
            this->orientation = RectanglePtr->orientation;
        }
        // else radius is left as it is
    }
    return *this;
};
Rectangle& Rectangle::operator=(const Rectangle& other) noexcept{
    if(this!= &other){
        Object::operator=(other);
            this->sides = other.sides;
            this->orientation = other.orientation;
    }
    return *this;
};  //uses Object copy operator= and copies radius
Rectangle& Rectangle::operator=(Rectangle&& other) noexcept{
    if(this!= &other){
        Object::operator=(std::move(other));
            this->sides = other.sides;
            this->orientation = other.orientation;
    }
    return *this;
};       //uses Object move operator= and copies radius

std::ostream& operator<<(std::ostream& out, const Rectangle& obj){ 
    out << static_cast<const Object&>(obj);
    out << " sides: " << obj.sides;
    out << " orientation: " << obj.orientation;
    return out;
};