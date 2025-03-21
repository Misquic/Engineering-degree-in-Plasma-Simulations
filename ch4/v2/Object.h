#ifndef OBJECT_H
#define OBJECT_H
#include <iostream>
#include <tuple>
#include <string>
#include "all.h"
#include "Vec3.h"

/*ideas:
    make movable and immovable objects
        -immovable: contributes as Dirichlet conditions
        -NOT DOING THIS ANYMORE movable: contributes as Dirichlet conditions but moves according to some pattern, e.g. constant gravity, these are big object
        that aren't moved due to E or B because they are too big, falling charged ball
    make different basic shapes using inheritance

    + operators works as intended they copy/move only shared attributes
*/
class Objects; //forward declaration for Object container named Objects

class Object{ //abstract class for other shapes
protected:
    std::string name = "Object";
    type_calc3 pos;    //position of center of Object
    //type_calc3 vel;    //velocity

    /*material properties*/
    type_calc phi = 0;


public:
    
    /*constructors*/
    //Object(type_calc3 pos, type_calc3 vel, type_calc phi); //normal constructor // needs to be changed
    Object(type_calc3 pos, type_calc phi); //normal constructor for immovable
    Object(const Object& other) noexcept;  //copying constructor
    Object(Object&& other) noexcept;       //moving constructor

    /*destructor*/
    virtual ~Object() noexcept = default;

    /*operators*/
    virtual Object& operator=(const Object& other) noexcept;  //copies pos and vel, doesn't change movable
    virtual Object& operator=(Object&& other) noexcept;       //moves pos and vel, doesn't change movable

    /*methods*/
    void setPhi(type_calc phi) noexcept;
    type_calc getPhi() noexcept;
    // virtual std::ostream print();
    virtual void print(std::ostream& out) const;
    
    virtual bool inObject(const type_calc3& x) const = 0;
    virtual void lineIntersect(const type_calc3& x1, const type_calc3& x2, type_calc& tp, type_calc3& pos, type_calc3& n) const = 0; //return tp, x_intersection, normal to surface at the point of intersection
    
    friend std::ostream& operator<<(std::ostream& out, const Object& obj);
    friend class Objects;
};


//////////////////////////////////////////////////////// SPHERE ////////////////////////////////////////////////////////

class Sphere: public Object{
protected:
    type_calc radius;
    type_calc r_squared;
public:
    /*constructors*/
    //Sphere(type_calc3 pos, type_calc3 vel, type_calc phi, type_calc radius);
    Sphere(type_calc3 pos, type_calc phi, type_calc radius);
    Sphere(const Sphere& other);
    Sphere(Sphere&& other);

    /*destructor*/
    virtual ~Sphere() noexcept = default;

    /*operators*/
    virtual Sphere& operator=(const Object& other) noexcept override;  //uses Object copy operator= and copies radius if argument is Sphere
    virtual Sphere& operator=(Object&& other) noexcept override;       //uses Object move operator= and copies radius if argument is Sphere
    virtual Sphere& operator=(const Sphere& other) noexcept;           //uses Object copy operator= and copies radius
    virtual Sphere& operator=(Sphere&& other) noexcept;                //uses Object move operator= and copies radius

    virtual void print(std::ostream& out) const override;
    virtual bool inObject(const type_calc3& x) const override;
    virtual void lineIntersect(const type_calc3& x1, const type_calc3& x2, type_calc& tp, type_calc3& pos, type_calc3& n) const override;


    friend std::ostream& operator<<(std::ostream& out, Sphere obj);
};


//////////////////////////////////////////////////////// RECTANGLE ////////////////////////////////////////////////////////

class Rectangle: public Object{
protected:
    type_calc3 sides; // side sizes
    type_calc3 half_sides;
    type_calc3 orientation;
    type_calc3 x_min, x_max; // position of minimum vertex, position of maximum vertex

    void find_n(int side, type_calc3& n) const; //used in line intersect 
public:
    /*constructors*/
    //Rectangle(type_calc3 pos, type_calc3 vel, type_calc phi, type_calc3 sides, type_calc3 orientation);
    Rectangle(type_calc3 pos, type_calc phi, type_calc3 sides, type_calc3 orientation); //TODO implement orientation
    Rectangle(type_calc3 pos, type_calc phi, type_calc3 sides);
    Rectangle(const Rectangle& other);
    Rectangle(Rectangle&& other);

    /*destructor*/
    virtual ~Rectangle() noexcept = default;

    /*operators*/
    virtual Rectangle& operator=(const Object& other) noexcept override;  //uses Object copy operator= and copies radius if argument is Rectangle
    virtual Rectangle& operator=(Object&& other) noexcept override;       //uses Object move operator= and copies radius if argument is Rectangle
    virtual Rectangle& operator=(const Rectangle& other) noexcept;        //uses Object copy operator= and copies radius
    virtual Rectangle& operator=(Rectangle&& other) noexcept;             //uses Object move operator= and copies radius

    virtual void print(std::ostream& out) const override;
    virtual bool inObject(const type_calc3& x) const override;
    virtual void lineIntersect(const type_calc3& x1, const type_calc3& x2, type_calc& t_entry, type_calc3& pos, type_calc3& n) const override; //return tp, point of intersection and normal to reflected surface

    friend std::ostream& operator<<(std::ostream& out, Rectangle obj);

};


#endif