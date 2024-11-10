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
void Sphere::lineIntersect(const type_calc3& x1, const type_calc3& x2, type_calc& t_entry, type_calc3& intersection_point, type_calc3& n) const{ 
    type_calc3 B = x2 - x1;
    type_calc3 A = x1 - this->pos;
    type_calc a = B*B;
    type_calc b = 2*(A*B);
    type_calc c = A*A - this->r_squared;
    type_calc det = b*b - 4*a*c;


    if( det < 0){
        t_entry = 0.5;
    }else{
        t_entry = (-b + std::sqrt(det))/(2*a);
        if (t_entry<0 || t_entry>1.0){
            t_entry = (-b - std::sqrt(det))/(2*a);
            if (t_entry<0 || t_entry>1.0){
                t_entry = 0.5;
                std::cerr << "Failed to find a line-sphere intersection!" << std::endl;
            }
        }
    }

    intersection_point = x1 + t_entry*B; // intersect pos on surface
    n = (intersection_point-this->pos).unit();
    

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
Rectangle::Rectangle(type_calc3 pos, type_calc phi, type_calc3 sides): Object(pos, phi), orientation{0,0,0}{
    this->sides = {fabs(sides[0]), fabs(sides[1]), fabs(sides[2])};
    half_sides = sides*0.5;
    name = "Rectangle";
    x_min = pos - half_sides;
    x_max = pos + half_sides;
};;

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
void Rectangle::lineIntersect(const type_calc3& x1, const type_calc3& x2, type_calc& t_entry, type_calc3& pos, type_calc3& n) const{ //only if x2 is inside the box
    type_calc3 A = x2 - x1; // if x2 == x1, A[i] = inf, min and < 0 covers it

    type_calc3 t_x_min = (x_min - x1)/A; //time of intersection with x_min sides
    type_calc3 t_x_max = (x_max - x1)/A; //time of intersection with x_max sides

    int possible_sides[3] = {0,1,2}; // x-, y-, z-, x+, y+, z+
    //std::cout << "min: " << t_x_min << " max: " <<  t_x_max << "\n";


///////////////////////////

    for(int i = 0; i < 3; i++){
        if(t_x_min[i] > t_x_max[i]){
            std::swap(t_x_min[i], t_x_max[i]);
            possible_sides[i] += 3; // to find side of intersection
        }
    }
    //std::cout << "min: " << t_x_min << " max: " <<  t_x_max << "\n";


    //find max and side:
    int sides[3] = {0,0,0}; 
    t_entry = t_x_min[0];
    sides[0] = possible_sides[0];
    int j = 0, k = 0, t = 0; //k number of intersecting sides if particle intersects with vertex or edge
    for(int i = 0; i < 2; i++){
        if(t_entry < t_x_min[i+1]){
            t_entry = t_x_min[i+1];
            j = i+1;
            k = 0;
        }
        else if( std::abs(t_entry - t_x_min[i+1]) < 1e-6){
            k++;
            t = i+1;

        }
    }
    //std::cout << t_x_min[2] - t_x_min[0] << "\n";

    find_n(possible_sides[j], n); // 0, 1, 2, 3, 4 or 5 -> (-1,0,0), (0,-1,0), (0,0,-1), (1,0,0), (0,1,0), (0,0,1)
    
    if(k != 0){
        if(k<2){
            find_n(possible_sides[t], n);
        }
        else{
            find_n(possible_sides[2], n);
            find_n(possible_sides[1], n);
        }
        n.normalise();
    }

    pos = x1 + t_entry*(x2 - x1);
};

void Rectangle::find_n(int side, type_calc3& n) const{
    switch(side){
        case 0:
            n += {-1,0,0};
            break;
        case 1:
            n += {0,-1,0};
            break;
        case 2:
            n += {0,0,-1};
            break;
        case 3:
            n += {1,0,0};
            break;
        case 4:
            n += {0,1,0};
            break;
        case 5:
            n += {0,0,1};
            break;
        default:
            std::cerr << "Rectangle::lineintersect default for side, we shouldn't be here\n";
            n = {0,0,0};
            break;
    }
}


