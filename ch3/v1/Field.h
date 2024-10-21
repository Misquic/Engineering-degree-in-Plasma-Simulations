#ifndef FIELD_H
#define FIELD_H
#include <vector>
#include <ostream>
#include <iostream>
#include <sstream>
#include "all.h"
#include "Vec3.h"
//using namespace std;

// #define FOR(i, n) for(int i = 0; i < this->n; i++) /* works but bad practise :(*/

template <class data_type>
class Field{
protected:
    std::vector<std::vector<std::vector<data_type>>> data; //3D container

public:
    const int ni, nj, nk;  //number of nodes in individual directions
    const int nn[3];       //number of nodes

    /*constructors*/
    //Field();                               //normal constructor
    Field(int ni, int nj, int nk);         //normal constructor
    Field(int nn[3]);         //normal constructor
    Field(int3 nn);         //normal constructor
    Field(const Field<data_type>& other);  //copying constructor
    Field(Field<data_type>&& other);       //moving constructor

    /*deconstructor*/
    ~Field(); //is it neccesary?

    /*methods*/
    std::string print() const;
    void clear();
    void scatter(type_calc3 lc, type_calc val);
    data_type gather(type_calc3 lc) const;

    /*operators []*/
    std::vector<std::vector<data_type>>& operator[](int i);
    const std::vector<std::vector<data_type>>& operator[](int i) const;
    //data_type operator()(int i, int j, int k) const; /*not needed if above is present?*/

    /*field-scalar operators*/
    Field<data_type>& operator=(const data_type val);
    Field<data_type> operator+(const data_type val) const;
    Field<data_type> operator-(const data_type val) const;
    Field<data_type> operator*(const data_type val) const;
    Field<data_type> operator/(const data_type val) const;
    void operator+=(const data_type val);
    void operator-=(const data_type val);
    void operator*=(const data_type val);
    void operator/=(const data_type val);

    /*field-field operators*/
    Field<data_type>& operator=(const Field<data_type>& other);  //copying operator
    Field<data_type>& operator=(Field<data_type>&& other);       //moving operator
    Field<data_type> operator+(const Field<data_type>& other) const;
    Field<data_type> operator-(const Field<data_type>& other) const;
    Field<data_type> operator*(const Field<data_type>& other) const;
    Field<data_type> operator/(const Field<data_type>& other) const;
    void operator+=(const Field<data_type>& other);
    void operator-=(const Field<data_type>& other);
    void operator*=(const Field<data_type>& other);
    void operator/=(const Field<data_type>& other);


};

/*constructors:*/
// template <class data_type>
// Field<data_type>::Field(): ni{0}, nj{0}, nk{0}, nn{ni, nj, nk} {
//     this->data = std::vector<std::vector<std::vector<data_type>>>(ni, std::vector<std::vector<data_type>>(nj, std::vector<data_type>(nk)));
// };
template <class data_type>
Field<data_type>::Field(int ni, int nj, int nk): ni{ni}, nj{nj}, nk{nk}, nn{ni, nj, nk} {
    this->data = std::vector<std::vector<std::vector<data_type>>>(ni, std::vector<std::vector<data_type>>(nj, std::vector<data_type>(nk)));
};
template <class data_type>
Field<data_type>::Field(int nn[3]): ni{nn[0]}, nj{nn[1]}, nk{nn[2]}, nn{ni, nj, nk}{
    this->data = std::vector<std::vector<std::vector<data_type>>>(ni, std::vector<std::vector<data_type>>(nj, std::vector<data_type>(nk)));
};
template <class data_type>
Field<data_type>::Field(int3 nn): ni{nn[0]}, nj{nn[1]}, nk{nn[2]}, nn{ni, nj, nk}{
    this->data = std::vector<std::vector<std::vector<data_type>>>(ni, std::vector<std::vector<data_type>>(nj, std::vector<data_type>(nk)));
};
template <class data_type>
Field<data_type>::Field(const Field<data_type>& other): ni{other.ni}, nj{other.nj}, nk{other.nk}, nn{other.ni, other.nj, other.nk} {
    //std::cout << "copied\n";
    this->data = std::vector<std::vector<std::vector<data_type>>>(ni, std::vector<std::vector<data_type>>(nj, std::vector<data_type>(nk)));
    
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] = other.data[i][j][k]; /*flag, maybe needed operator() as const?*//*no not needed*/
            }
        }
    }
};
template <class data_type>
Field<data_type>::Field(Field<data_type>&& other): ni{other.ni}, nj{other.nj}, nk{other.nk}, nn{other.ni, other.nj, other.nk} {
    std::cout << "moved\n";

    if(!(this->data.empty())) this->data.clear();
    this->data = std::vector<std::vector<std::vector<data_type>>>(ni, std::vector<std::vector<data_type>>(nj, std::vector<data_type>(nk)));
    this->data = move(other.data);
    other.data.clear();
};

/*destructor*/
template <class data_type>
Field<data_type>::~Field() {
    // Since std::vector handles memory management, we don't need to explicitly delete anything.
    // We just clear the data to ensure proper cleanup.
    if (!data.empty()) {
        data.clear();  // Clear the std::vector, which will deallocate the memory.
    }
}

/*methods*/
template <class data_type>
std::string Field<data_type>::print() const{
    std::stringstream out;
    //std::cout << field.nj << " :nj\n";
    //std::cout << field.nk << " :nk\n";
    for(int k=0; k<this->nk; k++){
        out << "\nk: " << k;
        //std::cout << "k: " << k << "\n";
        for(int j=0; j<this->nj; j++){
        out << "\n";
            //std::cout << "j: " << j << "\n";
            for(int i=0; i<this->ni; i++){
                //std::cout << "i: " << i << "\n";
                out << this->data[i][j][k] << " ";
            }
        }
    }
    return out.str();
};
template <class data_type>
void Field<data_type>::clear(){
    (*this) = 0;
};
template <class data_type> // in field because always used with fields
void Field<data_type>::scatter(type_calc3 lc, type_calc val){
    // if(lc[0] < 0 || lc[0] > ni-1 || lc[1] < 0 || lc[1] > nj-1 || lc[2] < 0 || lc[2] > nk-1){
    //     std::cerr << "scatter out of bands\n";
    //     return;
    // }
    int       i      = (int)lc[0];  //calculate cell indexes and fractional distances
    int       j      = (int)lc[1];
    int       k      = (int)lc[2];

    type_calc di     = lc[0] - i;   //precalculate values
    type_calc dj     = lc[1] - j;
    type_calc dk     = lc[2] - k;
    type_calc one_di = 1 - di;
    type_calc one_dj = 1 - dj;
    type_calc one_dk = 1 - dk;

    data[i][j][k]       += (data_type)val * (one_di) * (one_dj) * (one_dk);
    data[i][j][k+1]     += (data_type)val * (one_di) * (one_dj) * (dk);
    data[i][j+1][k]     += (data_type)val * (one_di) * (dj)     * (one_dk);
    data[i][j+1][k+1]   += (data_type)val * (one_di) * (dj)     * (dk);
    data[i+1][j][k]     += (data_type)val * (di)     * (one_dj) * (one_dk);
    data[i+1][j][k+1]   += (data_type)val * (di)     * (one_dj) * (dk);
    data[i+1][j+1][k]   += (data_type)val * (di)     * (dj)     * (one_dk);
    data[i+1][j+1][k+1] += (data_type)val * (di)     * (dj)     * (dk);
    
};
template <class data_type>
data_type Field<data_type>::gather(type_calc3 lc) const{
    data_type val{};

    int       i      = (int)lc[0];  //calculate cell indexes and fractional distances
    int       j      = (int)lc[1];
    int       k      = (int)lc[2];

    type_calc di     = lc[0] - i;   //precalculate values
    type_calc dj     = lc[1] - j;
    type_calc dk     = lc[2] - k;
    type_calc one_dj = 1 - dj;
    type_calc one_di = 1 - di;
    type_calc one_dk = 1 - dk;

    val = data[i][j][k] * one_di * one_dj * one_dk +
          data[i][j][k+1] * (one_di) * (one_dj) * (dk) + 
          data[i][j+1][k] * (one_di) * (dj) * (one_dk) + 
          data[i][j+1][k+1] * (one_di) * (dj) * (dk) + 
          data[i+1][j][k] * (di) * (one_dj) * (one_dk) + 
          data[i+1][j][k+1] * (di) * (one_dj) * (dk) + 
          data[i+1][j+1][k] * (di) * (dj) * (one_dk) + 
          data[i+1][j+1][k+1] * (di) * (dj) * (dk);
    return val;
};


/*operators[]*/
template <class data_type>
std::vector<std::vector<data_type>>& Field<data_type>::operator[](int i) { 
    //std::cout << "access\n";
    return data[i]; 
};
template <class data_type>
const std::vector<std::vector<data_type>>& Field<data_type>::operator[](int i) const{ 
    //std::cout << "const access\n";
    return data[i]; 
};

/*field-scalar operators*/
template <class data_type>
Field<data_type>& Field<data_type>::operator=(const data_type val){
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] = val;
            }
        }
    }
    return *this;
};
template <class data_type>
Field<data_type> Field<data_type>::operator+(const data_type val) const {
    Field<data_type> ret(this->ni, this->nj, this->nk);

    // int i,j,k;
    // FOR(i, ni){
    //     FOR(j, nj){
    //         FOR(k, nk){ /*works but bad practise :(*/
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] + val;
            }
        }
    }
    return ret;
};
template <class data_type>
Field<data_type> Field<data_type>::operator-(const data_type val) const {
    Field<data_type> ret(this->ni, this->nj, this->nk);

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] - val;
            }
        }
    }
    return ret;
};
template <class data_type>
Field<data_type> Field<data_type>::operator*(const data_type val) const {
    Field<data_type> ret(this->ni, this->nj, this->nk);

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] * val;
            }
        }
    }
    return ret;
};
template <class data_type>
Field<data_type> Field<data_type>::operator/(const data_type val) const {
    Field<data_type> ret(this->ni, this->nj, this->nk);
    const data_type inverse_val = 1.0/val;
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] * inverse_val;  /*is faster than dividing*/
            }
        }
    }
    return ret;
};
template <class data_type>
void Field<data_type>::operator+=(const data_type val){
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] += val;
            }
        }
    }
};
template <class data_type>
void Field<data_type>::operator-=(const data_type val){
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] -= val;
            }
        }
    }
};
template <class data_type>
void Field<data_type>::operator*=(const data_type val){
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] *= val;
            }
        }
    }
};
template <class data_type>
void Field<data_type>::operator/=(const data_type val){
    const data_type inverse_val = 1.0/val;
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] *= inverse_val; /*is faster than dividing*/
            }
        }
    }
};

/*field-field operators*/
template <class data_type>
Field<data_type>& Field<data_type>::operator=(const Field<data_type>& other){
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field = Field.");
        return *this;
    }
    
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] = other[i][j][k];
            }
        }
    }
    return *this;
};
template <class data_type>
Field<data_type>& Field<data_type>::operator=(Field<data_type>&& other){
    if(this != &other){
        if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
            throw std::invalid_argument("Dimensions of both fields must be the same for Field = &&Field.");
            return *this; /*seems like an arbitral choice?*/
        }

        if(!(this->data.empty())) this->data.clear();
        this->data = move(other.data);
        other.data.clear(); 
    }

    return *this;
};
template <class data_type>
Field<data_type> Field<data_type>::operator+(const Field<data_type>& other) const {
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field + Field.");
        return *this;
    }

    Field<data_type> ret(this->ni, this->nj, this->nk);
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] + other.data[i][j][k];
            }
        }
    }
    return ret;
};
template <class data_type>
Field<data_type> Field<data_type>::operator-(const Field<data_type>& other) const {
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field - Field.");
        return *this;
    }

    Field<data_type> ret(this->ni, this->nj, this->nk);
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] - other.data[i][j][k];
            }
        }
    }
    return ret;
};
template <class data_type>
Field<data_type> Field<data_type>::operator*(const Field<data_type>& other) const { /*most likely not needed operator*/
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field * Field.");
        return *this;
    }

    Field<data_type> ret(this->ni, this->nj, this->nk);
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] * other.data[i][j][k];
            }
        }
    }
    return ret;
};
template <class data_type>
Field<data_type> Field<data_type>::operator/(const Field<data_type>& other) const { /*most likely not needed operator*/
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field / Field.");
        return *this;
    }

    Field<data_type> ret(this->ni, this->nj, this->nk);
    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                ret.data[i][j][k] = this->data[i][j][k] / other.data[i][j][k];
            }
        }
    }
    return ret;
};
template <class data_type>
void Field<data_type>::operator+=(const Field<data_type>& other){
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field += Field.");
        return;
    }

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] += other[i][j][k];
            }
        }
    }
};
template <class data_type>
void Field<data_type>::operator-=(const Field<data_type>& other){
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field -= Field.");
        return;
    }

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] -= other[i][j][k];
            }
        }
    }
};
template <class data_type>
void Field<data_type>::operator*=(const Field<data_type>& other){
    if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field *= Field.");
        return;
    }

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                this->data[i][j][k] *= other[i][j][k];
            }
        }
    }
};
template <class data_type>
void Field<data_type>::operator/=(const Field<data_type>& other){
    if (ni != other.ni || nj != other.nj || nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field /= Field.");
        return;
    }

    for(int i = 0; i < ni; i++){
        for(int j = 0; j < nj; j++){
            for(int k = 0; k < nk; k++){
                if(other.data[i][j][k]!=0){
                    data[i][j][k] /= other.data[i][j][k];
                }
                else{
                    data[i][j][k] = 0;
                }
            }
        }
    }
};





/*sclalar-field operators*/
template <class data_type>
Field<data_type> operator+(const data_type& val, const Field<data_type>& field) {
    return field + val;
};
template <class data_type>
Field<data_type> operator-(const data_type& val, const Field<data_type>& field) {
    return field - val;
};
template <class data_type>
Field<data_type> operator*(const data_type& val, const Field<data_type>& field) {
    return field * val;
};
template <class data_type>
Field<data_type> operator/(const data_type& val, const Field<data_type>& field) {
    return field / val;
};

template<typename data_type, typename val_name>
Field<data_type> operator+(const val_name& val, const Field<data_type>& field) {
    return field + data_type( val);
};
template<typename data_type, typename val_name>
Field<data_type> operator-(const val_name& val, const Field<data_type>& field) {
    return field - data_type( val);
};
template<typename data_type, typename val_name>
Field<data_type> operator*(const val_name& val, const Field<data_type>& field) {
    return field * data_type( val);
};
template<typename data_type, typename val_name>
Field<data_type> operator/(const val_name& val, const Field<data_type>& field) {
    return field / data_type( val);
};

/*functions*/
template <class data_type>
std::ostream& operator<<(std::ostream &out, const Field<data_type> &field){
    //std::cout << field.nj << " :nj\n";
    //std::cout << field.nk << " :nk\n";
    for(int k=0; k<field.nk; k++){
        //out << "\nk: " << k;
        //std::cout << "k: " << k << "\n";
        for(int j=0; j<field.nj; j++){
            //out << "\n";
            //std::cout << "j: " << j << "\n";
            for(int i=0; i<field.ni; i++){
                //std::cout << "i: " << i << "\n";
                out << field[i][j][k] << " ";
            }
        }
        out << "\n";
    }
    return out;
};



#endif