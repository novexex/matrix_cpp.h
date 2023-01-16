#ifndef SRC_S21_MATRIX_OOP_H_
#define SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>
#include <cstring>
#include <exception>

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
    S21Matrix operator +(const S21Matrix &other)const;
    S21Matrix operator -(const S21Matrix &other)const;
    S21Matrix operator *(const double other)const;
    S21Matrix operator *(const S21Matrix &other)const;
    bool operator ==(const S21Matrix &other)const;
    double& operator() (int i, int j)const;

    ~S21Matrix();

    int get_rows()const;
    int get_columns()const;
    void set_rows(int i);
    void set_columns(int j);
    bool eq_matrix(const S21Matrix &other)const;
    void sum_matrix(const S21Matrix &other)const;
    void sub_matrix(const S21Matrix &other)const;
    void mul_number(const double num)const;
    void mul_matrix(const S21Matrix &other);
    S21Matrix transpose()const;
    S21Matrix calc_complements();
    double determinant();
    S21Matrix inverse_matrix();

 private:
    int rows_, columns_;
    double** matrix_;

    S21Matrix minor_matrix(int rows, int columns);
    void delete_matrix();
    void new_matrix();
    void set_matrix(int rows, int columns);
};

#endif  // SRC_S21_MATRIX_OOP_H_
