#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>
#include <cstring>
#include <exception>

#define s21_EPS 1e-7

class S21Matrix {
 public:
    S21Matrix();
    S21Matrix(int rows, int columns);
    S21Matrix(const S21Matrix &other);
    S21Matrix(S21Matrix &&other);

    S21Matrix& operator =(const S21Matrix &other);
    S21Matrix& operator +=(const S21Matrix &other);
    S21Matrix& operator -=(const S21Matrix &other);
    S21Matrix& operator *=(const S21Matrix &other);
    S21Matrix& operator *=(const double &num);
    S21Matrix operator +(const S21Matrix &other);
    S21Matrix operator -(const S21Matrix &other);
    S21Matrix operator *(const double other);
    S21Matrix operator *(const S21Matrix &other);
    bool operator ==(const S21Matrix &other);
    double& operator() (int i, int j);

    ~S21Matrix();

    int get_rows();
    int get_columns();
    void set_rows(int i);
    void set_columns(int j);
    bool eq_matrix(const S21Matrix &other);
    void sum_matrix(const S21Matrix &other);
    void sub_matrix(const S21Matrix &other);
    void mul_number(const double num);
    void mul_matrix(const S21Matrix &other);
    S21Matrix transpose();
    S21Matrix calc_complements();
    double determinant();
    S21Matrix inverse_matrix();

 private:
    int rows_;
    int columns_;
    double** matrix_;

    S21Matrix minor_matrix(int rows, int columns);
    void delete_matrix();
    void new_matrix();
    void set_matrix(int rows, int columns);
};

#endif  // SRC_S21_MATRIX_OOP_H_
