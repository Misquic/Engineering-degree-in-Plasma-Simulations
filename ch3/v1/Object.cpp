#include "Object.h"

/*constructors*/
Object::Object(type_calc3 pos, type_calc3 vel, bool isMovable, World& world): pos{pos}, vel{vel}, isMovable{isMovable}, world{world}{
    if(!isMovable){
        vel.clear();
    }
};
Object::Object(const Object& other) noexcept: pos{other.pos}, vel{other.vel}, isMovable{other.isMovable}, world{other.world}{ //does it work as intended?
};
Object::Object(Object&& other) noexcept: pos{std::move(other.pos)}, vel{std::move(other.vel)}, isMovable{other.isMovable}, world{other.world}{
};

/*destructors*/
// Object::~Object(){ //not needed here, because type_calc3 is type of Vec3 which has its own destructor
// }

/*operators*/
Object& Object::operator=(const Object& other) noexcept{
    if(this!= &other){
        this->pos = other.pos;
        this->vel = other.vel;
    }
    return *this;
};
Object& Object::operator=(Object&& other) noexcept{
    if(this!= &other){
        this->pos = std::move(other.pos);
        this->vel = std::move(other.vel);
    }
    return *this;
};

/*methods*/
void Object::advance(){
    if(isMovable){
        pos += vel * world.getDt(); //moves object
        
        //if forces to objects are coded
        //vel += acc / mass * world.getDt();
    }
};

std::ostream& operator<<(std::ostream& out, const Object& obj){
    out << "Name: " << std::setw(9) << obj.name;
    out << " pos: " << obj.pos;
    //if(obj.isMovable){
        out << " vel: " << obj.vel;
    //}
    return out;
};

//////////////////////////////////////////////////////// SPHERE ////////////////////////////////////////////////////////

/*constructors*/
Sphere::Sphere(type_calc3 pos, type_calc3 vel, bool isMovable, World& world, type_calc radius): Object(pos, vel, isMovable, world), radius{radius}{
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
Rectangle::Rectangle(type_calc3 pos, type_calc3 vel, bool isMovable, World& world, type_calc3 x0, type_calc3 x1, type_calc3 orientation):
     Object(pos, vel, isMovable, world), x0{x0}, x1{x1}, orientation{orientation}{
    name = "Rectangle";
};
Rectangle::Rectangle(const Rectangle& other): Object(other), x0{other.x0}, x1{other.x1}, orientation{other.orientation}{
    name = "Rectangle";
};
Rectangle::Rectangle(Rectangle&& other): Object(std::move(other)), x0{other.x0}, x1{other.x1}, orientation{other.orientation}{
    name = "Rectangle";
};

/*operators*/
Rectangle& Rectangle::operator=(const Object& other) noexcept{
    if(this!= &other){
        Object::operator=(other);
        if(const Rectangle* RectanglePtr = dynamic_cast<const Rectangle*>(&other)){
            this->x1 = RectanglePtr->x1;
            this->x0 = RectanglePtr->x0;
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
            this->x1 = RectanglePtr->x1;
            this->x0 = RectanglePtr->x0;
            this->orientation = RectanglePtr->orientation;
        }
        // else radius is left as it is
    }
    return *this;
};
Rectangle& Rectangle::operator=(const Rectangle& other) noexcept{
    if(this!= &other){
        Object::operator=(other);
            this->x1 = other.x1;
            this->x0 = other.x0;
            this->orientation = other.orientation;
    }
    return *this;
};  //uses Object copy operator= and copies radius
Rectangle& Rectangle::operator=(Rectangle&& other) noexcept{
    if(this!= &other){
        Object::operator=(std::move(other));
            this->x1 = other.x1;
            this->x0 = other.x0;
            this->orientation = other.orientation;
    }
    return *this;
};       //uses Object move operator= and copies radius

std::ostream& operator<<(std::ostream& out, const Rectangle& obj){ 
    out << static_cast<const Object&>(obj);
    out << " x0: " << obj.x0;
    out << " x1: " << obj.x1;
    out << " orientation: " << obj.orientation;
    return out;
};