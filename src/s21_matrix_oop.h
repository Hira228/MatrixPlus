#ifndef CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_

#include <algorithm>
#include <cmath>
#include <iostream>

#define EPS 1e-7

class S21Matrix {
 public:
  // Constructors and assignment operators
  S21Matrix() : S21Matrix(0, 0) {}
  explicit S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;

  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix operator=(S21Matrix&& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(double num);
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(const double num);
  friend S21Matrix operator*(const double num, const S21Matrix& other);
  bool operator==(const S21Matrix& other) const;
  double& operator()(int i, int j) noexcept;

  // Destructor
  ~S21Matrix() noexcept;

  // getters & setters
  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(int rows);
  void SetCols(int cols);

  // All other functions
  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

 private:
  int rows_ = 0, cols_ = 0;
  double** matrix_ = nullptr;

  // All other support functions
  void OutOfRange(int i, int j) const;
  void ExceptionMatrixDoesntExist() const;
  void ExceptionMatrixDoesntExist(const S21Matrix& other) const;
  void ExceptionMatrixIsNotSquare() const;
  void ExceptionDifferentDimensionsForMul(const S21Matrix& other) const;
  void ExceptionDifferentDimensions(const S21Matrix& other) const;
  void NewMatrixForMinor(S21Matrix& minor, int crossed_out_rows,
                         int crossed_out_column) noexcept;
  void AllocNewMatrix(double**& matrix, int rows, int cols) noexcept;
};

#endif  // CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_