#include "s21_matrix_oop.h"
#include <gtest/gtest.h>

TEST(Matrix, False) {
    S21Matrix matrix_x;
    S21Matrix matrix_a = std::move(matrix_x);
    S21Matrix matrix_b(3, 3);
    matrix_b(1, 1) = 2;
    ASSERT_FALSE(matrix_a == matrix_b);
}

TEST(Determinant, True) {
    S21Matrix matrix_a(3, 3);
    matrix_a(0, 0) = 2;
    matrix_a(0, 1) = 5;
    matrix_a(0, 2) = 7;
    matrix_a(1, 0) = 6;
    matrix_a(1, 1) = 3;
    matrix_a(1, 2) = 4;
    matrix_a(2, 0) = 5;
    matrix_a(2, 1) = -2;
    matrix_a(2, 2) = -3;
    ASSERT_EQ(matrix_a.determinant(),  -1);
}

TEST(Determinant1, True) {
    S21Matrix matrix_a(1, 1);
    matrix_a(0, 0) = 2;
    ASSERT_EQ(matrix_a.determinant(),  2);
}

TEST(Calc_complements, True) {
    S21Matrix matrix_a(3, 3);
    matrix_a(0, 0) = 1;
    matrix_a(0, 1) = 2;
    matrix_a(0, 2) = 3;
    matrix_a(1, 0) = 0;
    matrix_a(1, 1) = 4;
    matrix_a(1, 2) = 2;
    matrix_a(2, 0) = 5;
    matrix_a(2, 1) = 2;
    matrix_a(2, 2) = 1;

    S21Matrix matrix_b = matrix_a.calc_complements();

    S21Matrix matrix_c(3, 3);
    matrix_c(0, 0) = 0;
    matrix_c(0, 1) = 10;
    matrix_c(0, 2) = -20;
    matrix_c(1, 0) = 4;
    matrix_c(1, 1) = -14;
    matrix_c(1, 2) = 8;
    matrix_c(2, 0) = -8;
    matrix_c(2, 1) = -2;
    matrix_c(2, 2) = 4;
    ASSERT_TRUE(matrix_b == matrix_c);
}

TEST(Inverse_matrix, True) {
    S21Matrix matrix_a(3, 3);
    matrix_a(0, 0) = 2;
    matrix_a(0, 1) = 5;
    matrix_a(0, 2) = 7;
    matrix_a(1, 0) = 6;
    matrix_a(1, 1) = 3;
    matrix_a(1, 2) = 4;
    matrix_a(2, 0) = 5;
    matrix_a(2, 1) = -2;
    matrix_a(2, 2) = -3;

    S21Matrix matrix_b = matrix_a.inverse_matrix();

    S21Matrix matrix_c(3, 3);
    matrix_c(0, 0) = 1;
    matrix_c(0, 1) = -1;
    matrix_c(0, 2) = 1;
    matrix_c(1, 0) = -38;
    matrix_c(1, 1) = 41;
    matrix_c(1, 2) = -34;
    matrix_c(2, 0) = 27;
    matrix_c(2, 1) = -29;
    matrix_c(2, 2) = 24;
    ASSERT_TRUE(matrix_b == matrix_c);
}

TEST(EqMatrix, True) {
    S21Matrix matrix_a(3, 3);
    S21Matrix matrix_b(3, 3);
    ASSERT_TRUE(matrix_a == matrix_b);
}

TEST(EqMatrix, False) {
    S21Matrix matrix_a(3, 3);
    S21Matrix matrix_b(2, 2);
    ASSERT_FALSE(matrix_a == matrix_b);
}

TEST(SumMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3.14;
    matrix_a(0, 1) = 0.56;
    matrix_a(1, 0) = -69.3;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -78.14;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -0.3;
    matrix_b(1, 1) = 2;

    result(0, 0) = -75;
    result(0, 1) = 0.56;
    result(1, 0) = -69.6;
    result(1, 1) = 2;
    matrix_a.sum_matrix(matrix_b);
    ASSERT_TRUE(matrix_a == result);
}

TEST(SubMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -7;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -3;
    matrix_b(1, 1) = 2;

    result(0, 0) = 10;
    result(0, 1) = 2;
    result(1, 0) = -3;
    result(1, 1) = -2;

    matrix_a.sub_matrix(matrix_b);

    ASSERT_TRUE(matrix_a == result);
}

TEST(MulNumber, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = -78.14;
    matrix_a(0, 1) = 0;
    matrix_a(1, 0) = -0.3;
    matrix_a(1, 1) = 2;

    result(0, 0) = -781.4;
    result(0, 1) = 0;
    result(1, 0) = -3;
    result(1, 1) = 20;

    matrix_a.mul_number(10);

    ASSERT_TRUE(matrix_a == result);
}

TEST(MulMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6.6;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -7;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -3.5;
    matrix_b(1, 1) = 2;

    result(0, 0) = -28;
    result(0, 1) = 4;
    result(1, 0) = 46.2;
    result(1, 1) = 0;

    matrix_a.mul_matrix(matrix_b);

    ASSERT_TRUE(matrix_a == result);
}

TEST(Parentheses, True) {
    S21Matrix matrix_a(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6.6;
    matrix_a(1, 1) = 0;
    ASSERT_EQ(matrix_a(0, 1), 2);
}

TEST(Parentheses, False) {
    S21Matrix matrix_a(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6.6;
    matrix_a(1, 1) = 0;
    ASSERT_NE(matrix_a(0, 1), 10);
}

TEST(Plus, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3.14;
    matrix_a(0, 1) = 0.56;
    matrix_a(1, 0) = -69.3;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -78.14;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -0.3;
    matrix_b(1, 1) = 2;

    result(0, 0) = -75;
    result(0, 1) = 0.56;
    result(1, 0) = -69.6;
    result(1, 1) = 2;

    ASSERT_TRUE((matrix_a + matrix_b) == result);
}

TEST(Minus, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -7;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -3;
    matrix_b(1, 1) = 2;

    result(0, 0) = 10;
    result(0, 1) = 2;
    result(1, 0) = -3;
    result(1, 1) = -2;

    ASSERT_TRUE((matrix_a - matrix_b) == result);
}

TEST(OperatorMultiplyMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6.6;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -7;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -3.5;
    matrix_b(1, 1) = 2;

    result(0, 0) = -28;
    result(0, 1) = 4;
    result(1, 0) = 46.2;
    result(1, 1) = 0;

    ASSERT_TRUE((matrix_a * matrix_b) == result);
}

TEST(OperatorMultiplyNumber, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = -78.14;
    matrix_a(0, 1) = 0;
    matrix_a(1, 0) = -0.3;
    matrix_a(1, 1) = 2;

    result(0, 0) = -781.4;
    result(0, 1) = 0;
    result(1, 0) = -3;
    result(1, 1) = 20;

    ASSERT_TRUE((matrix_a * 10) == result);
}

TEST(OperatorMultiplyNum, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = -78.14;
    matrix_a(0, 1) = 0;
    matrix_a(1, 0) = -0.3;
    matrix_a(1, 1) = 2;

    result(0, 0) = -781.4;
    result(0, 1) = 0;
    result(1, 0) = -3;
    result(1, 1) = 20;

    ASSERT_TRUE((matrix_a * 10) == result);
}

TEST(OperatorEqMatrix, True) {
    S21Matrix matrix_b(1, 3);
    S21Matrix matrix_a = matrix_b;
    ASSERT_TRUE(matrix_a == matrix_b);
    }

    TEST(OperatorSumMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3.14;
    matrix_a(0, 1) = 0.56;
    matrix_a(1, 0) = -69.3;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -78.14;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -0.3;
    matrix_b(1, 1) = 2;

    result(0, 0) = -75;
    result(0, 1) = 0.56;
    result(1, 0) = -69.6;
    result(1, 1) = 2;

    ASSERT_TRUE((matrix_a += matrix_b) == result);
}

TEST(OperatorSubMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -7;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -3;
    matrix_b(1, 1) = 2;

    result(0, 0) = 10;
    result(0, 1) = 2;
    result(1, 0) = -3;
    result(1, 1) = -2;

    ASSERT_TRUE((matrix_a -= matrix_b) == result);
}

TEST(OperatorMulNumber, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = -78.14;
    matrix_a(0, 1) = 0;
    matrix_a(1, 0) = -0.3;
    matrix_a(1, 1) = 2;

    result(0, 0) = -781.4;
    result(0, 1) = 0;
    result(1, 0) = -3;
    result(1, 1) = 20;

    ASSERT_TRUE((matrix_a *= 10) == result);
}

TEST(OperatorMulMatrix, True) {
    S21Matrix matrix_a(2, 2);
    S21Matrix matrix_b(2, 2);
    S21Matrix result(2, 2);

    matrix_a(0, 0) = 3;
    matrix_a(0, 1) = 2;
    matrix_a(1, 0) = -6.6;
    matrix_a(1, 1) = 0;

    matrix_b(0, 0) = -7;
    matrix_b(0, 1) = 0;
    matrix_b(1, 0) = -3.5;
    matrix_b(1, 1) = 2;

    result(0, 0) = -28;
    result(0, 1) = 4;
    result(1, 0) = 46.2;
    result(1, 1) = 0;

    ASSERT_TRUE((matrix_a *= matrix_b) == result);
}

TEST(Transpose, True) {
    S21Matrix matrix_a(4, 3);
    S21Matrix result(3, 4);

    matrix_a(0, 0) = 7;
    matrix_a(0, 1) = -98;
    matrix_a(0, 2) = 0.5;
    matrix_a(1, 0) = 0;
    matrix_a(1, 1) = 5.4;
    matrix_a(1, 2) = 32;
    matrix_a(2, 0) = 3.12;
    matrix_a(2, 1) = 2323;
    matrix_a(2, 2) = 23;
    matrix_a(3, 0) = -78;
    matrix_a(3, 1) = 476.4;
    matrix_a(3, 2) = 21;

    result(0, 0) = 7;
    result(0, 1) = 0;
    result(0, 2) = 3.12;
    result(0, 3) = -78;
    result(1, 0) = -98;
    result(1, 1) = 5.4;
    result(1, 2) = 2323;
    result(1, 3) = 476.4;
    result(2, 0) = 0.5;
    result(2, 1) = 32;
    result(2, 2) = 23;
    result(2, 3) = 21;
    S21Matrix res = matrix_a.transpose();
    ASSERT_TRUE(res == result);
}

TEST(Get, True) {
    S21Matrix matrix_a(3, 3);

    matrix_a(0, 0) = 2;
    matrix_a(0, 1) = 5;
    matrix_a(0, 2) = 7;
    matrix_a(1, 0) = 6;
    matrix_a(1, 1) = 3;
    matrix_a(1, 2) = 4;
    matrix_a(2, 0) = 5;
    matrix_a(2, 1) = -2;
    matrix_a(2, 2) = -3;

    ASSERT_EQ(matrix_a.get_rows(), 3);
    ASSERT_EQ(matrix_a.get_columns(), 3);
}

TEST(Set, True) {
    S21Matrix matrix_a(3, 3);
    S21Matrix result(3, 2);

    matrix_a(0, 0) = 2;
    matrix_a(0, 1) = 5;
    matrix_a(0, 2) = 7;
    matrix_a(1, 0) = 6;
    matrix_a(1, 1) = 3;
    matrix_a(1, 2) = 4;
    matrix_a(2, 0) = 5;
    matrix_a(2, 1) = -2;
    matrix_a(2, 2) = -3;

    result(0, 0) = 2;
    result(0, 1) = 5;

    result(1, 0) = 6;
    result(1, 1) = 3;

    result(2, 0) = 5;
    result(2, 1) = -2;
    matrix_a.set_columns(2);

    ASSERT_TRUE(matrix_a == result);

    S21Matrix matrix_b(3, 3);
    S21Matrix result_b(2, 3);

    matrix_b(0, 0) = 2;
    matrix_b(0, 1) = 5;
    matrix_b(0, 2) = 7;
    matrix_b(1, 0) = 6;
    matrix_b(1, 1) = 3;
    matrix_b(1, 2) = 4;
    matrix_b(2, 0) = 5;
    matrix_b(2, 1) = -2;
    matrix_b(2, 2) = -3;

    result_b(0, 0) = 2;
    result_b(0, 1) = 5;
    result_b(0, 2) = 7;
    result_b(1, 0) = 6;
    result_b(1, 1) = 3;
    result_b(1, 2) = 4;

    matrix_b.set_rows(2);
    ASSERT_TRUE(matrix_b == result_b);
}

TEST(OperatorAssignment, Throw) {
    S21Matrix tmp(2, 2);
    EXPECT_THROW(tmp = tmp, std::invalid_argument);
}

TEST(OperatorSum, Throw) {
    S21Matrix tmp(2, 3);
    S21Matrix tmp1(4, 5);
    EXPECT_THROW(tmp+tmp1, std::invalid_argument);
}

TEST(OperatorSub, Throw) {
    S21Matrix tmp(2, 3);
    S21Matrix tmp1(4, 5);
    EXPECT_THROW(tmp-tmp1, std::invalid_argument);
}

TEST(OperatorMul, Throw) {
    S21Matrix tmp(2, 3);
    S21Matrix tmp1(4, 5);
    EXPECT_THROW(tmp-tmp1, std::invalid_argument);
}

TEST(OperatorIndex, Throw) {
    S21Matrix tmp(2, 3);
    EXPECT_THROW(tmp(4, 5), std::out_of_range);
}

TEST(MulMatrix, Throw) {
    S21Matrix tmp(2, 3);
    S21Matrix tmp1(4, 5);
    EXPECT_THROW(tmp.mul_matrix(tmp1), std::invalid_argument);
}

TEST(CalcComplements, Throw) {
    S21Matrix tmp(2, 3);
    EXPECT_THROW(tmp.calc_complements(), std::invalid_argument);
}

TEST(Determinant, Throw) {
    S21Matrix tmp(2, 3);
    EXPECT_THROW(tmp.determinant(), std::invalid_argument);
}

TEST(InverseMatrix, Throw) {
    S21Matrix tmp(1, 1);
    EXPECT_THROW(tmp.inverse_matrix(), std::invalid_argument);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
