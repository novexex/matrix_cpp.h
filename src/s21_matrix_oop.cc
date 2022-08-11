#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : S21Matrix(3, 3) {}

S21Matrix::S21Matrix(int rows, int columns) {
    if (!rows && !columns)
        throw std::invalid_argument("Cannot be assigned.");
    rows_ = rows;
    columns_ = columns;
    new_matrix();
}

S21Matrix::S21Matrix(const S21Matrix &other) : rows_(other.rows_), columns_(other.columns_) {
    matrix_ = new double*[rows_]();
    for (int i = 0; i < rows_; i++) {
        matrix_[i] = new double[columns_]();
        std::memcpy(matrix_[i], other.matrix_[i], columns_ * sizeof(double));
    }
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), columns_(other.columns_), matrix_(other.matrix_) {
    other.rows_ = 0;
    other.columns_ = 0;
    other.matrix_ = nullptr;
}

S21Matrix& S21Matrix::operator =(const S21Matrix &other) {
    if (this == &other)
        throw std::invalid_argument("Cannot be assigned.");
    delete_matrix();
    rows_ = other.rows_;
    columns_ = other.columns_;
    new_matrix();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++)
            matrix_[i][j] = other.matrix_[i][j];
    }
    return *this;
}

S21Matrix& S21Matrix::operator +=(const S21Matrix &other) {
    sum_matrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator -=(const S21Matrix &other) {
    sub_matrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator *=(const S21Matrix &other) {
    mul_matrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator *=(const double &num) {
    mul_number(num);
    return *this;
}

S21Matrix S21Matrix::operator +(const S21Matrix &other)const {
    S21Matrix result(*this);
    result.sum_matrix(other);
    return result;
}

S21Matrix S21Matrix::operator -(const S21Matrix &other)const {
    S21Matrix result(*this);
    result.sub_matrix(other);
    return result;
}

S21Matrix S21Matrix::operator *(const double other)const {
    S21Matrix result(*this);
    result.mul_number(other);
    return result;
}

S21Matrix S21Matrix::operator *(const S21Matrix &other)const {
    S21Matrix result(*this);
    result.mul_matrix(other);
    return result;
}

bool S21Matrix::operator ==(const S21Matrix &other)const {
    return eq_matrix(other);
}

double& S21Matrix::operator() (int i, int j)const {
    if (i < 0 || i >= rows_ || j < 0 || j >= columns_)
        throw std::out_of_range("Index outside the matrix.");
    return matrix_[i][j];
}

S21Matrix::~S21Matrix() {
    delete_matrix();
}

int S21Matrix::get_rows()const {
    return rows_;
}

int S21Matrix::get_columns()const {
    return columns_;
}

void S21Matrix::set_rows(int i) {
    if (!i)
        throw std::invalid_argument("Cannot be assigned.");
    set_matrix(i, columns_);
}

void S21Matrix::set_columns(int j) {
    if (!j)
        throw std::invalid_argument("Cannot be assigned.");
    set_matrix(rows_, j);
}

void S21Matrix::set_matrix(int rows, int columns) {
    S21Matrix tmp(*this);
    delete_matrix();
    rows_ = rows;
    columns_ = columns;
    new_matrix();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++)
            matrix_[i][j] = tmp.matrix_[i][j];
    }
}

bool S21Matrix::eq_matrix(const S21Matrix &other)const {
    bool return_value = true;
    if (rows_ == other.rows_ && columns_ == other.columns_) {
        for (int i = 0; i < rows_; i++) {
            for (int j = 0; j < columns_; j++) {
                if (fabs(matrix_[i][j] - other.matrix_[i][j]) > s21_EPS)
                    return_value = false;
            }
        }
    } else {
        return_value = false;
    }
    return return_value;
}

void S21Matrix::sum_matrix(const S21Matrix &other)const {
    if (rows_ != other.rows_ && columns_ != other.columns_)
        throw std::invalid_argument("Different dimensionality of matrices.");
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++)
            matrix_[i][j] += other.matrix_[i][j];
    }
}

void S21Matrix::sub_matrix(const S21Matrix &other)const {
    if (rows_ != other.rows_ && columns_ != other.columns_)
        throw std::invalid_argument("Different dimensionality of matrices.");
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++)
            matrix_[i][j] -= other.matrix_[i][j];
    }
}
void S21Matrix::mul_number(const double num)const {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++)
            matrix_[i][j] *= num;
    }
}

void S21Matrix::mul_matrix(const S21Matrix &other) {
    if (columns_ != other.rows_)
        throw std::invalid_argument
        ("The number of columns of the first matrix does not equal the number of rows of the second matrix.");
    S21Matrix result(rows_, other.columns_);
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < other.columns_; j++) {
            for (int k = 0; k < columns_; k++)
                result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
        }
    }
    *this = result;
}

S21Matrix S21Matrix::transpose()const {
    S21Matrix result(columns_, rows_);
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++)
            result.matrix_[j][i] = matrix_[i][j];
    }
    return result;
}

S21Matrix S21Matrix::calc_complements() {
    if (columns_ != rows_)
        throw std::invalid_argument("The matrix is not square.");
    S21Matrix result(rows_, columns_);
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++) {
            double determinant = 0;
            S21Matrix new_minor = minor_matrix(i, j);
            determinant = new_minor.determinant();
            determinant = determinant * pow(-1, i + j);
            result.matrix_[i][j] = determinant;
        }
    }
    return result;
}

S21Matrix S21Matrix::minor_matrix(int rows, int columns) {
    S21Matrix result(rows_-1, columns_-1);
    for (int i = 0, mini_i = 0; i < rows_; i++) {
        if (i == rows)
            continue;
        for (int j = 0, mini_j = 0; j < columns_; j++) {
            if (j == columns)
                continue;
            result.matrix_[mini_i][mini_j] = matrix_[i][j];
            mini_j++;
        }
        mini_i++;
    }
    return result;
}

double S21Matrix::determinant() {
    if (columns_ != rows_)
        throw std::invalid_argument("The matrix is not square.");
    double result = 0;
    if (rows_ == 1) {
        result = matrix_[0][0];
    } else if (rows_ == 2) {
        result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
    } else {
        double determinant = 0;
        for (int i = 0; i < rows_; i++) {
            S21Matrix new_minor = minor_matrix(0, i);
            result = new_minor.determinant();
            determinant += pow(-1, i) * matrix_[0][i] * result;
        }
        result = determinant;
    }
    return result;
}

S21Matrix S21Matrix::inverse_matrix() {
    double determinant = this->determinant();
    if (fabs(determinant) < s21_EPS)
        throw std::invalid_argument("Matrix determinant equals 0.");
    S21Matrix new_minor = calc_complements();
    S21Matrix transpose_matrix = new_minor.transpose();
    S21Matrix res = transpose_matrix * (1./determinant);
    return res;
}

void S21Matrix::delete_matrix() {
    if (matrix_ != nullptr) {
        for (int i = 0; i < rows_; i++)
            delete []matrix_[i];
        delete []matrix_;
    }
}

void S21Matrix::new_matrix() {
    matrix_ = new double*[rows_]();
    for (int i = 0; i < rows_; i++)
        matrix_[i] = new double[columns_]();
}
