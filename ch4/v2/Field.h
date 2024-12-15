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
    int avg_samples = 0;

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
    void scatter(type_calc3 lc, data_type val);
    data_type gather(type_calc3 lc) const;
    void updateAverage(const Field<data_type>& other);
    void updateMovingAverage(const Field<data_type>& other); //it's not a real moving average, just average from up to 21 samples but can have less samples.
    int U(int i, int j, int k) const;

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
    /*fild-other_type field operators*/
    template<class data_type2>
    friend Field<data_type2> operator/(const Field<data_type2>& left, const Field<type_calc>& right);
    // template<class data_type2>
    // Field<data_type> operator/(const Field<data_type2>& other) const;


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
    //std::cout << "moved\n";

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

    (*this) = (data_type){};
};
template <class data_type> // in field because always used with fields
void Field<data_type>::scatter(type_calc3 lc, data_type val){
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
    data_type v_one_di_one_dj = val * one_di * one_dj;
    data_type v_one_dj_dj = val * one_di * dj;
    data_type v_di_one_dj = val * di * one_dj;
    data_type v_di_dj = val * di * dj;
    std::vector<data_type>& data_ij = data[i][j];
    std::vector<data_type>& data_ij1 = data[i][j+1];
    std::vector<data_type>& data_i1j = data[i+1][j];
    std::vector<data_type>& data_i1j1 = data[i+1][j+1];

    // data[i][j][k]       += (data_type)val * (one_di) * (one_dj) * (one_dk);
    // data[i][j][k+1]     += (data_type)val * (one_di) * (one_dj) * (dk);
    // data[i][j+1][k]     += (data_type)val * (one_di) * (dj)     * (one_dk);
    // data[i][j+1][k+1]   += (data_type)val * (one_di) * (dj)     * (dk);
    // data[i+1][j][k]     += (data_type)val * (di)     * (one_dj) * (one_dk);
    // data[i+1][j][k+1]   += (data_type)val * (di)     * (one_dj) * (dk);
    // data[i+1][j+1][k]   += (data_type)val * (di)     * (dj)     * (one_dk);
    // data[i+1][j+1][k+1] += (data_type)val * (di)     * (dj)     * (dk);

    data_ij[k]       += v_one_di_one_dj * one_dk;
    data_ij[k+1]     += v_one_di_one_dj * dk;
    data_ij1[k]     += v_one_dj_dj     * one_dk;
    data_ij1[k+1]   += v_one_dj_dj     * dk;
    data_i1j[k]     += v_di_one_dj     * one_dk;
    data_i1j[k+1]   += v_di_one_dj     * dk;
    data_i1j1[k]   += v_di_dj         * one_dk;
    data_i1j1[k+1] += v_di_dj         * dk;
    
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
    type_calc one_di_one_dj = one_di * one_dj;
    type_calc one_dj_dj = one_di * dj;
    type_calc di_one_dj = di * one_dj;
    type_calc di_dj = di * dj;
    const std::vector<data_type>& data_ij = data[i][j];
    const std::vector<data_type>& data_ij1 = data[i][j+1];
    const std::vector<data_type>& data_i1j = data[i+1][j];
    const std::vector<data_type>& data_i1j1 = data[i+1][j+1];

    val = data_ij[k] * one_di_one_dj * one_dk +
          data_ij[k+1] * one_di_one_dj * dk + 
          data_ij1[k] * one_dj_dj * one_dk + 
          data_ij1[k+1] * one_dj_dj * dk + 
          data_i1j[k] * di_one_dj * one_dk + 
          data_i1j[k+1] * di_one_dj * dk + 
          data_i1j1[k] * di_dj * one_dk + 
          data_i1j1[k+1] * di_dj * dk;
    return val;
};
template <class data_type>
void Field<data_type>::updateAverage(const Field<data_type>& other){
    avg_samples++;

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                data[i][j][k] = (other[i][j][k] + data[i][j][k]*(avg_samples - 1.0))/avg_samples;
            }
        }
    }
};
template <class data_type>
void Field<data_type>::updateMovingAverage(const Field<data_type>& other){
    avg_samples++;
    if(avg_samples>20){
        clear();
        avg_samples = 1;
    };

    for(int i = 0; i < this->ni; i++){
        for(int j = 0; j < this->nj; j++){
            for(int k = 0; k < this->nk; k++){
                data[i][j][k] = (other[i][j][k] + data[i][j][k]*(avg_samples - 1.0))/avg_samples;
            }
        }
    }
    //data[0][0][0] = avg_samples;
};
template <class data_type>
int Field<data_type>::U(int i, int j, int k) const{
    return k*ni*nj + j*ni + i;
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
Field<data_type> Field<data_type>::operator*(const Field<data_type>& other) const {
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
                if(other.data[i][j][k]!=0){
                    ret.data[i][j][k] = this->data[i][j][k] / other.data[i][j][k];
                }
                else{
                    ret.data[i][j][k] = 0;
                }
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

/*field-other_type field operators*/
template <class data_type2>
Field<data_type2> operator/(const Field<data_type2>& left, const Field<type_calc>& other){
    if (left.ni != other.ni || left.nj != other.nj || left.nk != other.nk) {
        throw std::invalid_argument("Dimensions of both fields must be the same for Field<data_type> / Field<data_type2>.");
        return left;
    }

    Field<data_type2> ret(left.ni, left.nj, left.nk);
    for(int i = 0; i < left.ni; i++){
        for(int j = 0; j < left.nj; j++){
            for(int k = 0; k < left.nk; k++){
                
                if(other.data[i][j][k]!=0){
                    ret.data[i][j][k] = left.data[i][j][k] / other.data[i][j][k];
                }
                else{
                    ret.data[i][j][k] = 0;
                }
            }
        }
    }
    return ret;
};

// template<class data_type>
// template<class data_type2>
// Field<data_type> Field<data_type>::operator/(const Field<data_type2>& other) const{
//     if (this->ni != other.ni || this->nj != other.nj || this->nk != other.nk) {
//         throw std::invalid_argument("Dimensions of both fields must be the same for Field<data_type> / Field<data_type2>.");
//         return *this;
//     }

//     Field<data_type> ret(this->ni, this->nj, this->nk);
//     for(int i = 0; i < this->ni; i++){
//         for(int j = 0; j < this->nj; j++){
//             for(int k = 0; k < this->nk; k++){
//                 if(other.data[i][j][k]!=0){
//                     ret.data[i][j][k] = this->data[i][j][k] / other.data[i][j][k];
//                 }
//                 else{
//                     ret.data[i][j][k] = 0;
//                 }
//             }
//         }
//     }
//     return ret;
// };




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