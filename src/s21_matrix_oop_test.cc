#include "s21_matrix_oop.h"

#include <gtest/gtest.h>

TEST(Test, DefaultConstructor) {
  S21Matrix matrix;
  ASSERT_EQ(matrix.GetRows(), 0);
  ASSERT_EQ(matrix.GetCols(), 0);
}

TEST(Test, Constructor1) {
  S21Matrix matrix(2, 4);
  ASSERT_EQ(matrix.GetRows(), 2);
  ASSERT_EQ(matrix.GetCols(), 4);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 4; ++j) ASSERT_EQ(matrix(i, j), 0);
}

TEST(Test, Constructor2) {
  EXPECT_THROW(S21Matrix matrix(-2, 4), std::bad_alloc);
}

TEST(Test, CopyConstructor) {
  S21Matrix matrix(2, 3);
  matrix(0, 0) = 1;
  matrix(0, 1) = 2;
  matrix(0, 2) = 3.1;
  matrix(1, 0) = 4.5;
  matrix(1, 1) = 5;
  matrix(1, 2) = 6.2;
  S21Matrix copyMatrix = matrix;
  ASSERT_EQ(copyMatrix.GetRows(), 2);
  ASSERT_EQ(copyMatrix.GetCols(), 3);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) ASSERT_EQ(matrix(i, j), copyMatrix(i, j));
}

TEST(Test, MoveConstructor) {
  S21Matrix matrix(2, 3);
  matrix(0, 0) = 1;
  matrix(0, 1) = 2;
  matrix(0, 2) = 3.1;
  matrix(1, 0) = 4.5;
  matrix(1, 1) = 5;
  matrix(1, 2) = 6.2;
  S21Matrix matrixTemp(matrix);
  S21Matrix copyMatrix(std::move(matrix));
  ASSERT_EQ(copyMatrix.GetRows(), 2);
  ASSERT_EQ(copyMatrix.GetCols(), 3);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) ASSERT_EQ(matrixTemp(i, j), copyMatrix(i, j));
}

TEST(Test, EqMatrix1) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(1, 0) = 11.3;
  a(1, 1) = -1.3;
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 11.3;
  b(1, 1) = -1.3;
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(Test, EqMatrix2) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(1, 0) = 11.3;
  a(1, 1) = -1.3;
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 10;
  b(1, 1) = -1.3;
  EXPECT_FALSE(a.EqMatrix(b));
}

TEST(Test, EqMatrix3) {
  S21Matrix a(2, 3);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(0, 2) = 7;
  a(1, 0) = 11.3;
  a(1, 1) = -1.3;
  a(1, 2) = -1;
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 11.3;
  b(1, 1) = -1.3;
  EXPECT_FALSE(a.EqMatrix(b));
}

TEST(Test, SumMatrix1) {
  S21Matrix a(2, 3);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(0, 2) = 7;
  a(1, 0) = 11.3;
  a(1, 1) = -1.3;
  a(1, 2) = -1;
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 11.3;
  b(1, 1) = -1.3;
  EXPECT_THROW(a.SumMatrix(b), std::length_error);
}

TEST(Test, SumMatrix2) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(1, 0) = 1;
  a(1, 1) = -1;
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;
  a.SumMatrix(b);
  EXPECT_FALSE(a.EqMatrix(res));
}

TEST(Test, SumMatrix3) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(1, 0) = 0;
  a(1, 1) = -1;
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;
  a.SumMatrix(b);
  EXPECT_FALSE(a.EqMatrix(res));
}

TEST(Test, SubMatrix1) {
  S21Matrix a(2, 3);
  S21Matrix b(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(0, 2) = 7;
  a(1, 0) = 11.3;
  a(1, 1) = -1.3;
  a(1, 2) = -1;
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 11.3;
  b(1, 1) = -1.3;
  EXPECT_THROW(a.SumMatrix(b), std::length_error);
}

TEST(Test, SubMatrix2) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 1;
  a(1, 1) = -1;
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = -5;
  res(0, 1) = 3;
  res(1, 0) = 6;
  res(1, 1) = -7;
  a.SubMatrix(b);
  EXPECT_TRUE(a.EqMatrix(res));
}

TEST(Test, SubMatrix3) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 0;
  a(1, 1) = -1;
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = -5;
  res(0, 1) = 3;
  res(1, 0) = 6;
  res(1, 1) = -7;
  a.SubMatrix(b);
  EXPECT_FALSE(a.EqMatrix(res));
}

TEST(Test, MulNumber_1) {
  S21Matrix m(1, 1);
  m(0, 0) = 1;
  m = 5 * m;
  const double num = 5;
  S21Matrix res(1, 1);
  res(0, 0) = 5;
  EXPECT_TRUE(m.EqMatrix(res));
}

TEST(Test, MulNumber) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 0;
  a(1, 1) = -1;
  const double num = 2.5;
  S21Matrix res(2, 2);
  res(0, 0) = 10;
  res(0, 1) = 7.5;
  res(1, 0) = 0;
  res(1, 1) = -2.5;
  a.MulNumber(num);
  EXPECT_TRUE(a.EqMatrix(res));
}

TEST(Test, MulMatrix1) {
  S21Matrix a(3, 2);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  S21Matrix b(3, 4);
  b(0, 0) = 1;
  b(0, 1) = -1;
  b(0, 2) = 1;
  b(0, 3) = 4;
  b(1, 0) = 2;
  b(1, 1) = 3;
  b(1, 2) = 4;
  b(1, 3) = -5;
  b(2, 0) = 1;
  b(2, 1) = 1;
  b(2, 2) = 1;
  b(2, 3) = 1;

  EXPECT_THROW(a.MulMatrix(b), std::length_error);
}

TEST(Test, MulMatrix2) {
  S21Matrix a(3, 2);
  S21Matrix b(2, 3);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  b(0, 0) = 1;
  b(0, 1) = -1;
  b(0, 2) = 1;
  b(1, 0) = 2;
  b(1, 1) = 3;
  b(1, 2) = 4;
  S21Matrix res(3, 3);
  res(0, 0) = 9;
  res(0, 1) = 11;
  res(0, 2) = 17;
  res(1, 0) = 12;
  res(1, 1) = 13;
  res(1, 2) = 22;
  res(2, 0) = 15;
  res(2, 1) = 15;
  res(2, 2) = 27;
  a.MulMatrix(b);
  EXPECT_TRUE(a.EqMatrix(res));
}

TEST(Test, Transponse) {
  S21Matrix a(3, 2);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  S21Matrix res(2, 3);
  res(0, 0) = 1;
  res(0, 1) = 2;
  res(0, 2) = 3;
  res(1, 0) = 4;
  res(1, 1) = 5;
  res(1, 2) = 6;
  S21Matrix temp(2, 3);
  temp = a.Transpose();
  EXPECT_TRUE(res.EqMatrix(temp));
}

TEST(Test, CalcComplements1) {
  S21Matrix a(3, 3);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(1, 0) = 0;
  a(1, 1) = 4;
  a(1, 2) = 2;
  a(2, 0) = 5;
  a(2, 1) = 2;
  a(2, 2) = 1;
  S21Matrix res(3, 3);
  res(0, 0) = 0;
  res(0, 1) = 10;
  res(0, 2) = -20;
  res(1, 0) = 4;
  res(1, 1) = -14;
  res(1, 2) = 8;
  res(2, 0) = -8;
  res(2, 1) = -2;
  res(2, 2) = 4;
  EXPECT_TRUE(res.EqMatrix(a.CalcComplements()));
}

TEST(Test, Det1) {
  S21Matrix a(4, 3);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(1, 0) = 4;
  a(1, 1) = 5;
  a(1, 2) = 6;
  a(2, 0) = 7;
  a(2, 1) = 8;
  a(2, 2) = 9;
  a(3, 0) = 1;
  a(3, 1) = 1;
  a(3, 2) = 1;
  EXPECT_THROW(a.Determinant(), std::length_error);
}

TEST(Test, Det2) {
  S21Matrix a(3, 3);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(1, 0) = 4;
  a(1, 1) = 5;
  a(1, 2) = 6;
  a(2, 0) = 7;
  a(2, 1) = 8;
  a(2, 2) = 9;
  EXPECT_EQ(a.Determinant(), 0);
}

TEST(Test, Det3) {
  S21Matrix a(3, 3);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(1, 0) = 4;
  a(1, 1) = 5;
  a(1, 2) = 6;
  a(2, 0) = 7;
  a(2, 1) = 8;
  a(2, 2) = 10;
  EXPECT_DOUBLE_EQ(a.Determinant(), -3);
}

TEST(Test, Inverse1) {
  S21Matrix a(3, 3);
  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(1, 0) = 4;
  a(1, 1) = 5;
  a(1, 2) = 6;
  a(2, 0) = 7;
  a(2, 1) = 8;
  a(2, 2) = 9;
  EXPECT_THROW(a.InverseMatrix(), std::length_error);
}

TEST(Test, Inverse2) {
  S21Matrix a(3, 3);
  a(0, 0) = 2;
  a(0, 1) = 5;
  a(0, 2) = 7;
  a(1, 0) = 6;
  a(1, 1) = 3;
  a(1, 2) = 4;
  a(2, 0) = 5;
  a(2, 1) = -2;
  a(2, 2) = -3;
  S21Matrix res(3, 3);
  res(0, 0) = 1;
  res(0, 1) = -1;
  res(0, 2) = 1;
  res(1, 0) = -38;
  res(1, 1) = 41;
  res(1, 2) = -34;
  res(2, 0) = 27;
  res(2, 1) = -29;
  res(2, 2) = 24;
  EXPECT_TRUE(res.EqMatrix(a.InverseMatrix()));
}

TEST(Test, Eq1) {
  S21Matrix a(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(1, 0) = 11.3;
  a(1, 1) = -1.3;
  S21Matrix b(2, 2);
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 11.3;
  b(1, 1) = -1.3;
  EXPECT_TRUE(a == b);
}

TEST(Test, Eq2) {
  S21Matrix a(2, 2);
  a(0, 0) = 4.3;
  a(0, 1) = 3;
  a(1, 0) = 11.5;
  a(1, 1) = -1.3;
  S21Matrix b(2, 2);
  b(0, 0) = 4.3;
  b(0, 1) = 3;
  b(1, 0) = 11.3;
  b(1, 1) = -1.3;

  EXPECT_FALSE(a == b);
}

TEST(Test, Sum1) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 1;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;

  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;

  EXPECT_TRUE((a + b) == res);
}

TEST(Test, Sum2) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 9;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;

  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;

  EXPECT_FALSE((a + b) == res);
}

TEST(Test, Sum3) {
  S21Matrix a(1, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;

  EXPECT_THROW(a + b, std::length_error);
}

TEST(Test, Sub1) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 1;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;

  S21Matrix res(2, 2);
  res(0, 0) = -5;
  res(0, 1) = 3;
  res(1, 0) = 6;
  res(1, 1) = -7;

  EXPECT_TRUE((a - b) == res);
}

TEST(Test, Sub2) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 9;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;

  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;

  EXPECT_FALSE((a - b) == res);
}

TEST(Test, Sub3) {
  S21Matrix a(1, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;

  EXPECT_THROW(a - b, std::length_error);
}

TEST(Test, MulNumberOp) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 0;
  a(1, 1) = -1;
  const double num = 2.5;

  S21Matrix res(2, 2);
  res(0, 0) = 10;
  res(0, 1) = 7.5;
  res(1, 0) = 0;
  res(1, 1) = -2.5;

  EXPECT_TRUE(a * num == res);
}

TEST(Test, MulMatrixOp1) {
  S21Matrix a(3, 2);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  S21Matrix b(3, 4);
  b(0, 0) = 1;
  b(0, 1) = -1;
  b(0, 2) = 1;
  b(0, 3) = 4;
  b(1, 0) = 2;
  b(1, 1) = 3;
  b(1, 2) = 4;
  b(1, 3) = -5;
  b(2, 0) = 1;
  b(2, 1) = 1;
  b(2, 2) = 1;
  b(2, 3) = 1;
  EXPECT_THROW(a * b, std::length_error);
}

TEST(Test, MulMatrixOp2) {
  S21Matrix a(3, 2);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  S21Matrix b(2, 3);
  b(0, 0) = 1;
  b(0, 1) = -1;
  b(0, 2) = 1;
  b(1, 0) = 2;
  b(1, 1) = 3;
  b(1, 2) = 4;
  S21Matrix res(3, 3);
  res(0, 0) = 9;
  res(0, 1) = 11;
  res(0, 2) = 17;
  res(1, 0) = 12;
  res(1, 1) = 13;
  res(1, 2) = 22;
  res(2, 0) = 15;
  res(2, 1) = 15;
  res(2, 2) = 27;
  EXPECT_TRUE(a * b == res);
}

TEST(Test, LeftSum1) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 1;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;
  EXPECT_TRUE((a += b) == res);
}

TEST(Test, LeftSum3) {
  S21Matrix a(1, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  EXPECT_THROW(a += b, std::length_error);
}

TEST(Test, LeftSub1) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 1;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = -5;
  res(0, 1) = 3;
  res(1, 0) = 6;
  res(1, 1) = -7;
  EXPECT_TRUE((a -= b) == res);
}

TEST(Test, LeftSub2) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 9;
  a(1, 1) = -1;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  S21Matrix res(2, 2);
  res(0, 0) = 13;
  res(0, 1) = 3;
  res(1, 0) = -4;
  res(1, 1) = 5;
  EXPECT_FALSE((a -= b) == res);
}

TEST(Test, LeftSub3) {
  S21Matrix a(1, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  S21Matrix b(2, 2);
  b(0, 0) = 9;
  b(0, 1) = 0;
  b(1, 0) = -5;
  b(1, 1) = 6;
  EXPECT_THROW(a -= b, std::length_error);
}

TEST(Test, LeftMulNumberOp) {
  S21Matrix a(2, 2);
  a(0, 0) = 4;
  a(0, 1) = 3;
  a(1, 0) = 0;
  a(1, 1) = -1;
  const double num = 2.5;
  S21Matrix res(2, 2);
  res(0, 0) = 10;
  res(0, 1) = 7.5;
  res(1, 0) = 0;
  res(1, 1) = -2.5;
  EXPECT_TRUE((a *= num) == res);
}

TEST(Test, LeftMulMatrixOp1) {
  S21Matrix a(3, 2);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  S21Matrix b(3, 4);
  b(0, 0) = 1;
  b(0, 1) = -1;
  b(0, 2) = 1;
  b(0, 3) = 4;
  b(1, 0) = 2;
  b(1, 1) = 3;
  b(1, 2) = 4;
  b(1, 3) = -5;
  b(2, 0) = 1;
  b(2, 1) = 1;
  b(2, 2) = 1;
  b(2, 3) = 1;
  EXPECT_THROW(a *= b, std::length_error);
}

TEST(Test, LeftMulMatrixOp2) {
  S21Matrix a(3, 2);
  a(0, 0) = 1;
  a(0, 1) = 4;
  a(1, 0) = 2;
  a(1, 1) = 5;
  a(2, 0) = 3;
  a(2, 1) = 6;
  S21Matrix b(2, 3);
  b(0, 0) = 1;
  b(0, 1) = -1;
  b(0, 2) = 1;
  b(1, 0) = 2;
  b(1, 1) = 3;
  b(1, 2) = 4;
  S21Matrix res(3, 3);
  res(0, 0) = 9;
  res(0, 1) = 11;
  res(0, 2) = 17;
  res(1, 0) = 12;
  res(1, 1) = 13;
  res(1, 2) = 22;
  res(2, 0) = 15;
  res(2, 1) = 15;
  res(2, 2) = 27;
  EXPECT_TRUE((a *= b) == res);
}

TEST(Test, GetRows) {
  S21Matrix matrix(2, 2);
  EXPECT_TRUE(matrix.GetRows() == 2);
}

TEST(Test, GetCols) {
  S21Matrix matrix(2, 2);
  EXPECT_TRUE(matrix.GetRows() == 2);
}

TEST(Test, setRows) {
  S21Matrix matrix(2, 2);
  matrix.SetRows(5);
  EXPECT_TRUE(matrix.GetRows() == 5);
}

TEST(Test, setCols) {
  S21Matrix matrix(2, 2);
  matrix.SetCols(9);
  EXPECT_TRUE(matrix.GetCols() == 9);
}

TEST(Test, moveop) {
  S21Matrix a(2, 2);
  S21Matrix b;
  a(0, 0) = 4;
  a(0, 1) = 12;
  a(1, 0) = 1;
  a(1, 1) = 90;
  b = std::move(a);
  EXPECT_EQ(a.GetCols(), 0);
  EXPECT_EQ(a.GetRows(), 0);
  EXPECT_EQ(b(1, 1), 90);
}

TEST(Test, calccompserror) {
  S21Matrix mat(1, 2);
  EXPECT_THROW(mat.CalcComplements(), std::length_error);
}

TEST(Test, singulardeterminant) {
  S21Matrix mat(1, 1);
  mat(0, 0) = 5;
  EXPECT_DOUBLE_EQ(mat(0, 0), mat.Determinant());
}

TEST(Test, inverseerror) {
  S21Matrix mat(1, 2);
  EXPECT_THROW(mat.InverseMatrix(), std::length_error);
}

TEST(Test, selfassert) {
  S21Matrix mat(1, 1), mat2(1, 1);
  mat(0, 0) = 5;
  mat2(0, 0) = 5;
  mat = mat;
  EXPECT_DOUBLE_EQ(mat(0, 0), mat2(0, 0));
}

TEST(Test, copy1) {
  S21Matrix mat(5, 5), mat2(3, 3);
  mat = mat2;
  EXPECT_DOUBLE_EQ(mat(0, 0), mat2(0, 0));
}

TEST(Test, copy2) {
  S21Matrix mat(3, 3), mat2(5, 5);
  mat = mat2;
  EXPECT_DOUBLE_EQ(mat(0, 0), mat2(0, 0));
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}