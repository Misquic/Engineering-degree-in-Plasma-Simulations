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

class Object{
protected:
    std::string name = "Object";
    type_calc3 pos;    //position
    type_calc3 vel;    //velocity
    World&     world;  //reference to world object

public:
    const bool isMovable;
    
    /*constructors*/
    Object(type_calc3 pos, type_calc3 vel, bool isMovable, World& world); //normal constructor
    Object(const Object& other) noexcept;  //copying constructor, copies reference to world so that every object refers to the same world(hopefully)
    Object(Object&& other) noexcept;       //moving constructor, copies reference to world so that every object refers to the same world(hopefully)

    /*destructor*/
    virtual ~Object() noexcept = default;

    /*operators*/
    virtual Object& operator=(const Object& other) noexcept;  //copies pos and vel, doesn't change world and isMovable
    virtual Object& operator=(Object&& other) noexcept;       //moves pos and vel, doesn't change world and isMovable

    /*methods*/
    void advance();

    friend std::ostream& operator<<(std::ostream& out, const Object& obj);

};


//////////////////////////////////////////////////////// SPHERE ////////////////////////////////////////////////////////

class Sphere: public Object{
protected:
    type_calc radius;
public:
    /*constructors*/
    Sphere(type_calc3 pos, type_calc3 vel, bool isMovable, World& world, type_calc radius);
    Sphere(const Sphere& other);
    Sphere(Sphere&& other);

    /*destructor*/
    virtual ~Sphere() noexcept = default;

    /*operators*/
    virtual Sphere& operator=(const Object& other) noexcept override;  //uses Object copy operator= and copies radius if argument is Sphere
    virtual Sphere& operator=(Object&& other) noexcept override;       //uses Object move operator= and copies radius if argument is Sphere
    virtual Sphere& operator=(const Sphere& other) noexcept;           //uses Object copy operator= and copies radius
    virtual Sphere& operator=(Sphere&& other) noexcept;                //uses Object move operator= and copies radius

    friend std::ostream& operator<<(std::ostream& out, const Sphere& obj);

};


//////////////////////////////////////////////////////// RECTANGLE ////////////////////////////////////////////////////////

class Rectangle: public Object{
protected:
    type_calc3 x0;
    type_calc3 x1;
    type_calc3 orientation;
public:
    /*constructors*/
    Rectangle(type_calc3 pos, type_calc3 vel, bool isMovable, World& world, type_calc3 x0, type_calc3 x1, type_calc3 orientation = {1, 1, 1});
    Rectangle(const Rectangle& other);
    Rectangle(Rectangle&& other);

    /*destructor*/
    virtual ~Rectangle() noexcept = default;

    /*operators*/
    virtual Rectangle& operator=(const Object& other) noexcept override;  //uses Object copy operator= and copies radius if argument is Rectangle
    virtual Rectangle& operator=(Object&& other) noexcept override;       //uses Object move operator= and copies radius if argument is Rectangle
    virtual Rectangle& operator=(const Rectangle& other) noexcept;        //uses Object copy operator= and copies radius
    virtual Rectangle& operator=(Rectangle&& other) noexcept;             //uses Object move operator= and copies radius

    friend std::ostream& operator<<(std::ostream& out, const Rectangle& obj);
};


#endif