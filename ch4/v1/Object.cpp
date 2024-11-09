#include "Object.h"
#include <cmath>
#include <algorithm>

/*constructors*/
Object::Object(type_calc3 pos, type_calc3 vel, type_calc phi): pos{pos}, vel{vel}, phi{phi}, movable{true}{
};
Object::Object(type_calc3 pos, type_calc phi): pos{pos}, vel{vel}, phi{phi}, movable{false}{
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
type_calc Object::getPhi() noexcept{
    return phi;
};

bool Object::isMovable() noexcept{
    return movable;
}
void Object::print(std::ostream& out) const{
    out << "Name: " << std::setw(9) << name;
    out << " pos: " << pos;
    if(movable){
        out << " vel: " << vel;
    }
};

/*friends*/
std::ostream& operator<<(std::ostream& out, const Object& obj){
    obj.print(out);
    return out;
};

//////////////////////////////////////////////////////// SPHERE ////////////////////////////////////////////////////////

/*constructors*/
Sphere::Sphere(type_calc3 pos, type_calc3 vel, type_calc phi, type_calc radius): Object(pos, vel, phi), radius{fabs(radius)}{
    r_squared = radius*radius;
    name = "Sphere";
};
Sphere::Sphere(type_calc3 pos, type_calc phi, type_calc radius): Object(pos, phi), radius{fabs(radius)}{
    r_squared = radius*radius;
    name = "Sphere";
};
Sphere::Sphere(const Sphere& other): Object(other), radius{other.radius}{
    r_squared = radius*radius;
    name = "Sphere";
};
Sphere::Sphere(Sphere&& other): Object(std::move(other)), radius{other.radius}{
    r_squared = radius*radius;
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

void Sphere::print(std::ostream& out) const{ // for 
    Object::print(out);
    out << " radius: " << radius;
};
bool Sphere::inObject(type_calc3 x) const{
    type_calc3 r = x - pos;
    if(r * r <= r_squared) return true;
    return false;
};
type_calc Sphere::lineIntersect(const type_calc3& x1, const type_calc3& x2) const{
    type_calc3 B = x2 - x1;
    type_calc3 A = x1 - this->pos;
    type_calc a = B*B;
    type_calc b = 2*(A*B);
    type_calc c = A*A - this->r_squared;
    type_calc det = b*b - 4*a*c;

    if( det < 0) return 0.5;
    type_calc tp = (-b + std::sqrt(det))/(2*a);
    if (tp<0 || tp>1.0){
        tp = (-b - std::sqrt(det))/(2*a);
        if (tp<0 || tp>1.0){
            tp = 0.5;
            std::cerr << "Failed to find a line-sphere intersection!" << std::endl;
        }
    }
    return tp;

};


std::ostream& operator<<(std::ostream& out, Sphere obj){ // for std::cout << <Sphere>
    obj.print(out);
    return out;
};


//////////////////////////////////////////////////////// Rectangle ////////////////////////////////////////////////////////

/*constructors*/
Rectangle::Rectangle(type_calc3 pos, type_calc3 vel, type_calc phi, type_calc3 sides, type_calc3 orientation): 
 Object(pos, vel, phi), orientation{orientation}{
    this->sides = {fabs(sides[0]), fabs(sides[1]), fabs(sides[2])};
    half_sides = sides*0.5;
    name = "Rectangle";
};
Rectangle::Rectangle(type_calc3 pos, type_calc phi, type_calc3 sides, type_calc3 orientation): 
 Object(pos, phi), orientation{orientation}{
    this->sides = {fabs(sides[0]), fabs(sides[1]), fabs(sides[2])};
    half_sides = sides*0.5;
    name = "Rectangle";
    x_min = pos - half_sides;
    x_max = pos + half_sides;
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

void Rectangle::print(std::ostream& out) const{
    Object::print(out);
    out << " sides: " << sides;
    out << " orientation: " << orientation;
};

std::ostream& operator<<(std::ostream& out, Rectangle obj){
    obj.print(out);
    return out;
};
bool Rectangle::inObject(type_calc3 x) const{ //for now orientation is not taken into account, to do that convert to inertial system of Rectangle and than use normal isObject?
    //if(orientation == type_calc3(1, 0, 0))
    type_calc3 temp = abs(x - pos);
    for(int i = 0; i < 3; i++){
        if(temp[i] > half_sides[i]) return false;
    }
    return true;
};
type_calc Rectangle::lineIntersect(const type_calc3& x1, const type_calc3& x2) const{ //only if x2 is inside the box
    type_calc3 A = x2 - x1; // if x2 == x1, A[i] = inf, min and < 0 covers it

    type_calc3 t_x_min = (x_min - x1)/A; //time of intersection with x_min sides
    type_calc3 t_x_max = (x_max - x1)/A; //time of intersection with x_max sides

    //std::cout << "min: " << t_x_min << " max: " <<  t_x_max << "\n";

///////////////////////////

    for(int i = 0; i < 3; i++){
        if(t_x_min[i] > t_x_max[i]){
            std::swap(t_x_min[i], t_x_max[i]);
        }
    }

    type_calc t_entry = std::max({t_x_min[0], t_x_min[1], t_x_min[2]});
    return t_entry;
    type_calc t_exit = std::min({t_x_max[0], t_x_max[1], t_x_max[2]});

///////////////////////////



    // checks which plane is encountered first, but it doenst have to be side :(
    // type_calc t_mins[3]; //times of intersections with real sides that line marked 
    // for(int i = 0; i < 3; i++){
    //     if(t_x_min[i]<0){ //if time < 0 it means that x_min is "behind" particle so other side must be intersected, always at least 3 sides are in front of particle 
    //         t_mins[i] = t_x_max[i];
    //     }
    //     else{ //if in particular i-th direction two sides are in front of particle line
    //         t_mins[i] = std::min(t_x_min[i], t_x_max[i]);
    //     }
    // }

    // return std::min({t_mins[0], t_mins[1], t_mins[2]}); //we have now granted that times are for nearest sides in front of particle, we chose that one which is nearest

};



