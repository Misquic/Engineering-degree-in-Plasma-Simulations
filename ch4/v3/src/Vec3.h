#ifndef VEC3_H
#define VEC3_H

#include <ostream>
#include <iostream>
#include <typeinfo>
#include <iomanip>
#include <sstream>
#include <cmath>
#include "all.h"
//using namespace std;

//TODO add casting from int to double fo example

template <class data_type>
class Vec3{
protected:
    data_type m_data[3];
public:
    /*constructors*/
    Vec3() noexcept;                                                         //normal constructor
    Vec3(const data_type x, const data_type y, const data_type z) noexcept;  //normal constructor
    Vec3(const data_type other_data[3]) noexcept;                            //normal constructor
    Vec3(const data_type val) noexcept;                                      //normal constructor
    Vec3(const Vec3<data_type>& other) noexcept;                             //copying constructor
    Vec3(Vec3<data_type>&& other) noexcept;                                  //moving constructor

    /*deconstructor*/
    //~Vec3();

    /*methods*/
    void clear();
    type_calc length() const; //mag
    void normalise(); //
    Vec3<data_type> unit() const;
    Vec3<data_type> cross(const Vec3<data_type>& other) const;
    void rotate(Vec3<data_type> axis_unit_vec, type_calc sin, type_calc cos);
    Vec3<data_type> elWiseMult(const Vec3<data_type>& other) const;  //elementwise multiplication
    data_type volume() const;
    bool isNan() const;
    bool isNan(const std::string& message) const;

    /*operators[]*/
    data_type& operator[](int i);
    const data_type& operator[](int i) const;

    /*Vec3-scalar operators*/
    Vec3<data_type>& operator=(const data_type& val) noexcept;
    Vec3<data_type> operator+(const data_type val) const;
    Vec3<data_type> operator-(const data_type val) const;
    Vec3<data_type> operator*(const data_type val) const;
    Vec3<data_type> operator/(const data_type val) const;
    void operator+=(const data_type val);
    void operator-=(const data_type val);
    void operator*=(const data_type val);
    void operator/=(const data_type val);

    /*Vec3-Vec3 operators*/
    Vec3<data_type>& operator=(const Vec3<data_type>& other) noexcept;  //copying operator
    Vec3<data_type>& operator=(Vec3<data_type>&& other) noexcept;       //moving operator
    Vec3<data_type> operator+(const Vec3<data_type>& other) const;
    Vec3<data_type> operator-(const Vec3<data_type>& other) const;
    data_type       operator*(const Vec3<data_type>& other) const;  //dot multiplication
    Vec3<data_type> operator/(const Vec3<data_type>& other) const;  //elementwise division
    void operator+=(const Vec3<data_type>& other);
    void operator-=(const Vec3<data_type>& other);
    void operator*=(const Vec3<data_type>& other) = delete;
    void operator/=(const Vec3<data_type>& other) = delete;

    bool operator==(const Vec3<data_type>& other) const;

};

/*constructors*/
template <class data_type>
Vec3<data_type>::Vec3() noexcept: m_data{0,0,0}{
};
template <class data_type>
Vec3<data_type>::Vec3(const data_type x, const data_type y, const data_type z) noexcept: m_data{x, y, z}{
};
template <class data_type>
Vec3<data_type>::Vec3(const data_type other_data[3]) noexcept: m_data{other_data[0], other_data[1], other_data[2]}{
};
template <class data_type>
Vec3<data_type>::Vec3(const data_type val) noexcept: m_data{val, val, val}{
};
template <class data_type>
Vec3<data_type>::Vec3(const Vec3<data_type>& other) noexcept: m_data{other.m_data[0], other.m_data[1], other.m_data[2]}{
};
template <class data_type>
Vec3<data_type>::Vec3(Vec3<data_type>&& other) noexcept: m_data{other.m_data[0], other.m_data[1], other.m_data[2]}{
};

/*methods*/
template <class data_type>
void Vec3<data_type>::clear(){
    for(int i = 0; i < 3; i ++){
        m_data[i] = 0;
    }
};
template <class data_type>
type_calc Vec3<data_type>::length() const{
    return std::sqrt((*this)*(*this));
};
template <class data_type>
void Vec3<data_type>::normalise(){
    (*this)/=(*this).length();
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::unit() const{
    type_calc l = length();
    #ifdef DEBUG
    if(!l){
        std::cerr << "Vec3<" << typeid(data_type).name()<< ">::unit() division by zero\n";
    }
    #endif
    return Vec3<data_type>( (*this)/l );

    //return Vec3<data_type>( {0,0,0} ); 
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::cross(const Vec3<data_type>& other) const{
    return Vec3<data_type>( m_data[1]*other.m_data[2] - m_data[2]*other.m_data[1] , m_data[2]*other.m_data[0] - m_data[0]*other.m_data[2] , m_data[0]*other.m_data[1] - m_data[1]*other.m_data[0] );
};

template <class data_type>
void Vec3<data_type>::rotate(Vec3<data_type> axis_unit_vec, type_calc sin, type_calc cos){
    Vec3<data_type> v_perrallel = ((*this)*axis_unit_vec)*axis_unit_vec;
    Vec3<data_type> v_perpendicular = (*this)-v_perrallel;
    v_perpendicular = cos*v_perpendicular + sin * (axis_unit_vec.cross(*this));
    (*this) = v_perrallel + v_perpendicular;
};


/*operators[]*/
template <class data_type>
data_type& Vec3<data_type>::operator[](int i){
    return m_data[i];
};
template <class data_type>
const data_type& Vec3<data_type>::operator[](int i) const {
    return m_data[i];
}

/*vec3-scalar operators*/
template <class data_type>
Vec3<data_type>& Vec3<data_type>::operator=(const data_type& val) noexcept{
    for(int i = 0; i < 3; i++){
        m_data[i] = val;
    }
    return (*this); //nie było tego i było okej?
    
};

template <class data_type>
Vec3<data_type> Vec3<data_type>::operator+(const data_type val) const{
    Vec3<data_type> ret(*this); /*probably faster than creating empty vec3 and then adding, other behaviour for field (i think so?)*/
    for(int i = 0; i < 3; i++){
        ret.m_data[i] += val;
    }
    return ret;
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::operator-(const data_type val) const{
    Vec3<data_type> ret(*this);
    for(int i = 0; i < 3; i++){
        ret.m_data[i] -= val;
    }
    return ret;
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::operator*(const data_type val) const{
    Vec3<data_type> ret(*this);
    for(int i = 0; i < 3; i++){
        ret.m_data[i] *= val;
    }
    return ret;
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::operator/(const data_type val) const{
    Vec3<data_type> ret(*this);
    #ifdef DEBUG
    if(!val){
        std::cerr << "Vec3<" << typeid(data_type).name() << ">::operator/(const data_type val) Division by zero\n";
    }
    #endif
    data_type inverse_val = 1.0/val;
    for(int i = 0; i < 3; i++){
        ret.m_data[i] *= inverse_val;
    }
    return ret;
};
template <class data_type>
void Vec3<data_type>::operator+=(const data_type val){
    for(int i = 0; i < 3; i++){
        m_data[i] += val;
    }
};
template <class data_type>
void Vec3<data_type>::operator-=(const data_type val){
    for(int i = 0; i < 3; i++){
        m_data[i] -= val;
    }
};
template <class data_type>
void Vec3<data_type>::operator*=(const data_type val){
    for(int i = 0; i < 3; i++){
        m_data[i] *= val;
    }
};
template <class data_type>
void Vec3<data_type>::operator/=(const data_type val){
    #ifdef DEBUG
    if(!val){
        std::cerr << "Vec3<" << typeid(data_type).name() << ">::operator/=(const data_type val) Division by zero\n";
    }
    #endif
    data_type inverse_val = 1.0/val;
    for(int i = 0; i < 3; i++){
        m_data[i] *= inverse_val;
    }
};

/*vec3-vec3 operators*/
template <class data_type>
Vec3<data_type>& Vec3<data_type>::operator=(const Vec3<data_type>& other) noexcept{
    for(int i = 0; i < 3; i++){
        this->m_data[i] = other.m_data[i];
    }
    return *this;
};
template <class data_type>
Vec3<data_type>& Vec3<data_type>::operator=(Vec3<data_type>&& other) noexcept{
    if(this != &other){
        m_data[0] = std::move(other.m_data[0]);
        m_data[1] = std::move(other.m_data[1]);
        m_data[2] = std::move(other.m_data[2]);
        
    }
    return *this;
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::operator+(const Vec3<data_type>& other) const{
    Vec3<data_type> ret(*this);
    for(int i = 0; i < 3; i++){
        ret.m_data[i] += other.m_data[i];
    }
    ///
    //Vec3<data_type> ret();
    //for(int i = 0; i<3; i++){
    //    ret.m_data[i] = this->m_data[i] + other.m_data[i] /*which is faster?*/
    //}
    //
    ///
    return ret;
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::operator-(const Vec3<data_type>& other) const{
    Vec3<data_type> ret(*this);
    for(int i = 0; i < 3; i++){
        ret.m_data[i] -= other.m_data[i];
    }
    return ret;
};
template <class data_type>
data_type Vec3<data_type>::operator*(const Vec3<data_type>& other) const{ //point multiplication
    return m_data[0]*other.m_data[0] + m_data[1]*other.m_data[1] + m_data[2]*other.m_data[2] ;
};
template <class data_type>
Vec3<data_type> Vec3<data_type>::operator/(const Vec3<data_type>& other) const{ //elementwise division
    Vec3<data_type> ret(*this);
    for(int i = 0; i < 3; i++){
        #ifdef DEBUG
        if(!other.m_data[i]){
            std::cerr << "Vec3<" << typeid(data_type).name() << ">::operator/(const Vec3<data_type>& other) Division by zero\n";
        }
        #endif
        ret.m_data[i] /= other.m_data[i];
    }
    return ret;
}; 
template<class data_type>
Vec3<data_type> Vec3<data_type>::elWiseMult(const Vec3<data_type>& other) const{  //elementwise multiplication
    return {m_data[0]*other.m_data[0], m_data[1]*other.m_data[1], m_data[2]*other.m_data[2]};
}; 
template<class data_type>
data_type Vec3<data_type>::volume() const{
    return m_data[0]*m_data[1]*m_data[2];
};
template<class data_type>
bool Vec3<data_type>::isNan() const{
    for(int i = 0; i < 3; i++){
        if(std::isnan(m_data[i])){
            std::cerr << "vec3: " << *this << " nan i: " << i << " ";
            return true;
        }
    }
    return false;
};
template<class data_type>
bool Vec3<data_type>::isNan(const std::string& message) const{
    for(int i = 0; i < 3; i++){
        if(std::isnan(m_data[i])){
            std::cerr << message << " vec3: " << *this << " nan i: " << i << " ";
            return true;
        }
    }
    return false;
};


template <class data_type>
void Vec3<data_type>::operator+=(const Vec3<data_type>& other){
    for(int i = 0; i < 3; i++){
        this->m_data[i] += other.m_data[i];
    }
};
template <class data_type>
void Vec3<data_type>::operator-=(const Vec3<data_type>& other){
    for(int i = 0; i < 3; i++){
        this->m_data[i] -= other.m_data[i];
    }
};

template <class data_type>
bool Vec3<data_type>::operator==(const Vec3<data_type>& other) const{
    for(int i = 0; i < 3; i++){
        if(this->m_data[i] != other.m_data[i]) return false;
    }
    return true;
};

/*sclalar-vec3 operators*/
template <class data_type>
Vec3<data_type> operator+(const data_type& val, const Vec3<data_type>& vec3) {
    return vec3 + val;
};
template <class data_type>
Vec3<data_type> operator-(const data_type& val, const Vec3<data_type>& vec3) {
    return vec3 - val;
};
template <class data_type>
Vec3<data_type> operator*(const data_type& val, const Vec3<data_type>& vec3) {
    return vec3 * val;
};
// template <class data_type>
// Vec3<data_type> operator/(const data_type& val, const Vec3<data_type>& vec3) {
//      //not usede ever, done it on mind autopilot
//     return vec3 / val;
// };

template<typename data_type, typename val_name>
Vec3<data_type> operator+(const val_name& val, const Vec3<data_type>& vec3) {
    return vec3 + data_type( val);
};
template<typename data_type, typename val_name>
Vec3<data_type> operator-(const val_name& val, const Vec3<data_type>& vec3) {
    return vec3 - data_type( val);
};
template<typename data_type, typename val_name>
Vec3<data_type> operator*(const val_name& val, const Vec3<data_type>& vec3) {
    return vec3 * data_type( val);
};
// template<typename data_type, typename val_name>
// Vec3<data_type> operator/(const val_name& val, const Vec3<data_type>& vec3) {
//     //not needed ever, done it on mind autopilot
//    return vec3 / data_type( val);
// };

/*functions*/
template <class data_type>
std::ostream& operator<<(std::ostream& out, const Vec3<data_type>& vec3){
    out << std::setw(5) << vec3[0] << " " << std::setw(5) << vec3[1] << " " << std::setw(5) << vec3[2];
    return out;
};
template <class data_type>
Vec3<data_type> cross(const Vec3<data_type>& l, const Vec3<data_type>& r){
    return Vec3<data_type>( l[1]*r[2] - l[2]*r[1] , l[2]*r[0] - l[0]*r[2] , l[0]*r[1] - l[1]*r[0] );
};

template <class data_type>
Vec3<data_type>  abs(const Vec3<data_type> vec){
    return Vec3<data_type>(std::fabs(vec[0]), std::fabs(vec[1]), std::fabs(vec[2]));
};


using double3 = Vec3<double>;
using type_calc3 = Vec3<type_calc>;
using int3 = Vec3<int>;

#endif