#ifndef MATRIX_H
#define MATRIX_H
#include "all.h"
#include <vector>
#include <array>
#include <memory>

using tcvector = std::vector<type_calc>;
//TODO in v2 change approach to this?, why *rows instead of 

template <int S>
class Row{
public:
    type_calc a[S]; // coefficients
    int col[S]; // collumn indexes, or -1 if not set

    /*constructor*/
    Row() noexcept;
    Row(const Row<S>& other) noexcept;
    Row(Row<S>&& other) noexcept;

    /*operators*/
    // void    operator=(const Row<S>& other) noexcept;  //copying operator
    // void    operator=(Row<S>&& other) noexcept;       //moving operator
    Row<S>& operator=(const Row<S>& other) noexcept;  //copying operator
    Row<S>& operator=(Row<S>&& other) noexcept;       //moving operator

};

class Matrix{
public:
    static constexpr int nvals = 7; //number of non zero columns per row;
    const int n_unknowns; //number of unknown values
    
    /*constructors*/
    Matrix(int n_unknowns) noexcept;
    Matrix(const Matrix& other) noexcept;
    Matrix(Matrix&& other) noexcept;

    /*operators*/
    Matrix& operator=(const Matrix& other);
    Matrix& operator=(Matrix&& other);
    tcvector operator*(const tcvector& v) const; // matrix * vector = vector
    type_calc& operator()(int r, int c); // reference to A[r,c]
    const type_calc& operator()(int r, int c) const;
    type_calc get(int r, int c) const; // value at A[r,c]
    
    /*methods*/
    void clearRow(int r);
    Matrix diagSubtract(const tcvector& v);
    Matrix invDiagonal();
    type_calc multiplyRow(int r, const tcvector& v) const;
    //tcvector getRow(int r);    
    

protected:
    std::vector<Row<nvals>> rows;
};

///////////////////////////////// TCVECTOR ///////////////////////////////////////////
tcvector operator-(const tcvector& left, const tcvector& right) noexcept;
tcvector operator+(const tcvector& left, const tcvector& right) noexcept;
type_calc operator*(const tcvector& left, const tcvector& right) noexcept;

tcvector operator*(const type_calc val, const tcvector& v) noexcept;
tcvector operator*(const tcvector& v, const type_calc val) noexcept;

#endif