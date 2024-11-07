#include "Matrix.h"
#include <iostream>
#include <cassert>

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
    // Matrix ret(*this); //copy
    // for(int i = 0; i < n_unknowns; i++){
    //     ret(i,i) = (*this)(i,i) - v[i];
    // }
    // return ret;
    Matrix M(*this);	//make a copy
	for (int u=0;u<n_unknowns;u++) M(u,u)=(*this)(u,u)-v[u];
	return M;
};
Matrix Matrix::invDiagonal() const{
    Matrix ret(this->n_unknowns);
    for(int i = 0; i < n_unknowns; i++){
        ret(i,i) = 1/(*this)(i,i);
    }
    return ret;
};
type_calc Matrix::multiplyRow(int r, const tcvector& v) const{
    const Row<nvals>& row = rows[r];
    type_calc sum = 0;
    for(int i = 0; i < nvals; i++){
        if(row.col[i]>=0) sum+=row.a[i]*v[row.col[i]];
        else break;
    }
    return sum;
};


///////////////////////////////// TCVECTOR ///////////////////////////////////////////
tcvector operator-(const tcvector& left, const tcvector& right) noexcept{
    assert(left.size()==right.size());
    tcvector ret(left.size());
    for(int i = 0; i < left.size(); i++){
        ret[i] = left[i] - right[i];
    }
    return ret;
};
tcvector operator+(const tcvector& left, const tcvector& right) noexcept{
    // if(left.size()!=right.size()){
    //     throw std::invalid_argument("tcvectors must have the same dimensions, err4");
    // }
    assert(left.size()==right.size());
    tcvector ret(left.size());
    for(int i = 0; i < left.size(); i++){
        ret[i] = left[i] + right[i];
    }
    return ret;
};
type_calc operator*(const tcvector& left, const tcvector& right) noexcept{
    // if(left.size()!=right.size()){
    //     throw std::invalid_argument("tcvectors must have the same dimensions, err4");
    // }
    assert(left.size()==right.size());

    type_calc ret{};
    for(int i = 0; i < left.size(); i++){
        ret+= left[i] * right[i];
    }
    return ret;
};
tcvector operator*(const type_calc val, const tcvector& v) noexcept{
    tcvector ret(v.size());
    for(int i = 0; i < v.size(); i++){
        ret[i] = v[i] * val;
    }
    return ret;
};
tcvector operator*(const tcvector& v, const type_calc val) noexcept{
    return val * v;
};