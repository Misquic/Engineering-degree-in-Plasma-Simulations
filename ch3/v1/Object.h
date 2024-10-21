#ifndef OBJECT_H
#define OBJECT_H
#include <ostream>
#include <string>
#include "all.h"
#include "Vec3.h"
#include "World.h"

/*ideas:
    make movable and immovable objects
        -immovable: contributes as Dirichlet conditions
        -movable: contributes as Dirichlet conditions but moves according to some pattern, e.g. constant gravity, these are big object
        that aren't moved due to E or B because they are too big, falling charged ball
    make different basic shapes using inheritance

    + operators works as intended they copy/move only shared attributes
    TODO sugarcubing, adding objects to world taking them into account in Solver
*/
class Objects; //forward declaration for Object container named Objects

class Object{
protected:
    std::string name = "Object";
    type_calc3 pos;    //position
    type_calc3 vel;    //velocity

    /*material properties*/
    type_calc phi = 0;
    bool movable;

public:
    
    /*constructors*/
    Object(type_calc3 pos, type_calc3 vel); //normal constructor // needs to be changed
    Object(type_calc3 pos); //normal constructor for immovable
    Object(const Object& other) noexcept;  //copying constructor
    Object(Object&& other) noexcept;       //moving constructor

    /*destructor*/
    virtual ~Object() noexcept = default;

    /*operators*/
    virtual Object& operator=(const Object& other) noexcept;  //copies pos and vel, doesn't change movable
    virtual Object& operator=(Object&& other) noexcept;       //moves pos and vel, doesn't change movable

    /*methods*/
    void setPhi(type_calc phi) noexcept;
    bool isMovable() noexcept;
    // virtual std::ostream print();
    virtual void print(std::ostream& out) const;

    friend std::ostream& operator<<(std::ostream& out, const Object& obj);
    friend class Objects;
};


//////////////////////////////////////////////////////// SPHERE ////////////////////////////////////////////////////////

class Sphere: public Object{
protected:
    type_calc radius;
public:
    /*constructors*/
    Sphere(type_calc3 pos, type_calc3 vel, type_calc radius);
    Sphere(type_calc3 pos, type_calc radius);
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

    friend std::ostream& operator<<(std::ostream& out, Sphere obj);
};


//////////////////////////////////////////////////////// RECTANGLE ////////////////////////////////////////////////////////

class Rectangle: public Object{
protected:
    type_calc3 sides; // side sizes
    type_calc3 orientation;
public:
    /*constructors*/
    Rectangle(type_calc3 pos, type_calc3 vel, type_calc3 sides, type_calc3 orientation);
    Rectangle(type_calc3 pos, type_calc3 sides, type_calc3 orientation);
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

    friend std::ostream& operator<<(std::ostream& out, Rectangle obj);

};


#endif