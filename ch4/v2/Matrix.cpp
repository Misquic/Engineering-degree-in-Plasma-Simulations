#include "Matrix.h"
#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>

////////////////////////////////// ROW /////////////////////////////////////////////

template <int S>
Row<S>::Row() noexcept{
    for(int i = 0; i < S; i++){
        a[i] = 0;
        col[i] = -1;
    }
};
template <int S>
Row<S>::Row(const Row<S>& other) noexcept{
    for(int i = 0; i < S; i++){
        a[i] = other.a[i];
        col[i] = other.col[i];
    }
};
template <int S>
Row<S>::Row(Row<S>&& other) noexcept{
    for (int i = 0; i < S; ++i) {
        a[i] = std::move(other.a[i]);
        col[i] = std::move(other.col[i]);
    }
};

// template <int S>
// void Row<S>::operator=(const Row<S>& other) noexcept{
//     for(int i = 0; i < S; i++){
//         a[i] = other.a[i];
//         col[i] = other.col[i];
//     }
// };
// template <int S>
// void Row<S>::operator=(Row<S>&& other) noexcept{
//     for(int i = 0; i < S; i++){
//         a = std::move(other.a);
//         col = std::move(other.col);
//     }
// };
template <int S>
Row<S>& Row<S>::operator=(const Row<S>& other) noexcept{
    if (this != &other) {
        for(int i = 0; i < S; i++){
            a[i] = other.a[i];
            col[i] = other.col[i];
        }
    }
    return (*this);
};
template <int S>
Row<S>& Row<S>::operator=(Row<S>&& other) noexcept{
    if (this != &other) {
        for(int i = 0; i < S; i++){
            a[i] = std::move(other.a[i]);
            col[i] = std::move(other.col[i]);
            other.col[i] = -1;
        }
    }
    return (*this);
};

///////////////////////////////// MATRIX ///////////////////////////////////////////

/*constructors*/
Matrix::Matrix(int n_unknowns) noexcept: n_unknowns{n_unknowns}, rows{std::vector<Row<nvals>>(n_unknowns)}{
};
Matrix::Matrix(const Matrix& other) noexcept: Matrix(other.n_unknowns){
    rows = other.rows;
};
Matrix::Matrix(Matrix&& other) noexcept: Matrix(other.n_unknowns){
    rows = std::move(other.rows);
};
/*operators*/
Matrix& Matrix::operator=(const Matrix& other){
    if(this!= &other){
        if(this->n_unknowns!= other.n_unknowns){
            throw std::invalid_argument("number of unknowns in matrixes dont match, err1");
        }
        this->rows = other.rows;
    }
    return (*this);
};
Matrix& Matrix::operator=(Matrix&& other){
    if(this!= &other){
        if(this->n_unknowns!= other.n_unknowns){
            throw std::invalid_argument("number of unknowns in matrixes dont match, err2");
        }
        this->rows = std::move(other.rows);
    }
    return (*this);
};
tcvector Matrix::operator*(const tcvector& v) const{
    tcvector ret(n_unknowns);
    if(this->n_unknowns!= v.size()){
        throw std::invalid_argument("number of unknowns in matrixes dont match dimension of vector, err3");
    }
    for(int u = 0; u < n_unknowns; u ++){
        const Row<nvals>& row = rows[u]; //change to pointer?
        ret[u] = 0;
        for(int i = 0; i < nvals; i++){
            if(row.col[i] >=0 ){
                if(std::isnan(row.a[i])|| std::isnan(v[row.col[i]])){
                    throw std::runtime_error("matrix operator* nan\n");
                }
                
                ret[u]+= row.a[i] * v[row.col[i]];
            }
            else break; //end at the first -1 in col
        }
    }
    return ret;
};
type_calc& Matrix::operator()(int r, int c){
    Row<nvals>& row = rows[r];
    int i = 0;
    for(i = 0; i < nvals; i++){
        if(row.col[i] == c) break;
        if(row.col[i] < 0){
            row.col[i] = c;
            break;
        }
    }
    assert(i!=nvals);
    return row.a[i];
};
const type_calc& Matrix::operator()(int r, int c) const{
    const Row<nvals>& row = rows[r];
    int i = 0;
    for(i = 0; i < nvals; i++){
        if(row.col[i] == c) break;
    }
    assert(i!=nvals);
    return row.a[i];
};

/*methods*/
void Matrix::clearRow(int r){
    rows[r] = Row<nvals>();
};
Matrix Matrix::diagSubtract(const tcvector& v) const{
    Matrix M(*this);	//make a copy
	for(size_t i = 0; i < n_unknowns; i++){
        #ifdef DEBUG
        if(std::isnan((*this)(i,i) - v[i])){
            throw std::runtime_error("Matrix diagSubtract nan");
        }
        #endif
        M(i,i) = (*this)(i,i) - v[i];
    }
    return M;
};
Matrix Matrix::invDiagonal() const{
    Matrix ret(this->n_unknowns);
    for(size_t i = 0; i < n_unknowns; i++){
        ret(i,i) = 1/(*this)(i,i);
    }
    return ret;
};
type_calc Matrix::multiplyRow(int r, const tcvector& v) const{
    const Row<nvals>& row = rows[r];
    type_calc sum = 0;
    for(size_t i = 0; i < nvals; i++){
        if(row.col[i]>=0) sum+=row.a[i]*v[row.col[i]];
        else break;
    }
    return sum;
};


///////////////////////////////// TCVECTOR ///////////////////////////////////////////
tcvector operator-(const tcvector& left, const tcvector& right) noexcept{
    size_t n_unknowns = left.size();
    assert(n_unknowns==right.size());
    tcvector ret(n_unknowns);
    for(size_t i = 0; i < n_unknowns; i++){
        ret[i] = left[i] - right[i];
    }
    return ret;
};
tcvector operator+(const tcvector& left, const tcvector& right) noexcept{
    size_t n_unknowns = left.size();
    assert(n_unknowns==right.size());
    tcvector ret(n_unknowns);
    for(size_t i = 0; i < n_unknowns; i++){
        ret[i] = left[i] + right[i];
    }
    return ret;
};
type_calc operator*(const tcvector& left, const tcvector& right){
    size_t n_unknowns = left.size();
    assert(n_unknowns==right.size());
    type_calc ret{};
    for(size_t i = 0; i < n_unknowns; i++){
        #ifdef DEBUG
        if(std::isnan(left[i]) || std::isnan(right[i])){
            std::stringstream sstr;
            sstr << "type_calc operator*(const tcvector& left, const tcvector& right) left: " << left[i] << " right: " << right[i] << " i: " << i << "\n";
            throw std::runtime_error(sstr.str());
        }
        #endif
        ret += left[i] * right[i];
    }
    return ret;
};

tcvector operator*(const type_calc val, const tcvector& v) noexcept{
    size_t n_unknowns = v.size();
    tcvector ret(n_unknowns);
    for(size_t i = 0; i < n_unknowns; i++){
        ret[i] = v[i] * val;
    }
    return ret;
};
tcvector operator*(const tcvector& v, const type_calc val) noexcept{
    return val * v;
};