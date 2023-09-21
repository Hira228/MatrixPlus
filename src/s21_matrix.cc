#include "s21_matrix_oop.h"

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  try {
    matrix_ = new double* [rows] {};
  } catch (...) {
    delete[] matrix_;
    throw;
  }
  for (int i = 0; i < rows; ++i) {
    try {
      matrix_[i] = new double[cols]{};
    } catch (...) {
      for (int j = 0; j < i; ++j) {
        delete[] matrix_[j];
      }
      delete[] matrix_;
      throw;
    }
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : S21Matrix(other.rows_, other.cols_) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept {
  if (this != &other) {
    for (int i = 0; i < rows_; ++i) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }
}

S21Matrix::~S21Matrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  S21Matrix copy_matrix(other);
  *this = std::move(copy_matrix);
  return *this;
}

S21Matrix S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    for (int i = 0; i < rows_; ++i) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix copy_matrix = *this;
  copy_matrix.SumMatrix(other);
  return copy_matrix;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix copy_matrix = *this;
  copy_matrix.SubMatrix(other);
  return copy_matrix;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix copy_matrix = *this;
  copy_matrix.MulMatrix(other);
  return copy_matrix;
}

S21Matrix S21Matrix::operator*(double num) {
  S21Matrix copy_matrix = *this;
  copy_matrix.MulNumber(num);
  return copy_matrix;
}

S21Matrix operator*(const double num, const S21Matrix& other) {
  S21Matrix copy_matrix = other;
  copy_matrix.MulNumber(num);
  return copy_matrix;
}

S21Matrix& S21Matrix::operator*=(double num) {
  MulNumber(num);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(other);
}

double& S21Matrix::operator()(int i, int j) noexcept {
  OutOfRange(i, j);
  return matrix_[i][j];
}

int S21Matrix::GetRows() const noexcept { return rows_; }

int S21Matrix::GetCols() const noexcept { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows == rows_) return;
  if (rows < 1) {
    throw std::length_error("out of range");
  }
  double** newMatrix;
  AllocNewMatrix(newMatrix, rows, cols_);
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;

  rows_ = rows;
  matrix_ = newMatrix;
}

void S21Matrix::SetCols(int cols) {
  if (cols == cols_) return;
  if (cols < 1) {
    throw std::length_error("out of range");
  }

  double** newMatrix;
  AllocNewMatrix(newMatrix, rows_, cols);
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;

  cols_ = cols;
  matrix_ = newMatrix;
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  bool result = true;
  ExceptionMatrixDoesntExist(other);
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    result = false;
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        result = false;
        break;
      }
    }
  }
  return result;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  ExceptionMatrixDoesntExist(other);
  ExceptionDifferentDimensions(other);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  ExceptionMatrixDoesntExist(other);
  ExceptionDifferentDimensions(other);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  ExceptionMatrixDoesntExist();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  ExceptionMatrixDoesntExist(other);
  ExceptionDifferentDimensionsForMul(other);
  if (other.rows_ != cols_) {
    throw;
  }
  S21Matrix result_matrix(rows_, other.cols_);
  for (int i = 0; i < result_matrix.rows_; ++i) {
    for (int j = 0; j < result_matrix.cols_; ++j) {
      for (int k = 0; k < cols_; ++k) {
        result_matrix.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = result_matrix;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  ExceptionMatrixDoesntExist();
  ExceptionMatrixIsNotSquare();
  S21Matrix result(rows_, cols_);
  S21Matrix minor(rows_ - 1, cols_ - 1);
  double result_determinant = 0;

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      double alg_add = 1.0f;
      NewMatrixForMinor(minor, i, j);
      if ((i + j) % 2 == 1) {
        alg_add = -1.0f;
      }
      result_determinant = minor.Determinant();
      if (result_determinant == 0.0f) alg_add = 1;
      result.matrix_[i][j] = alg_add * result_determinant;
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  ExceptionMatrixIsNotSquare();
  int alg_add = 1;
  double result = 0.0f;
  if (rows_ == 1) {
    result = matrix_[0][0];
  } else if (rows_ == 2) {
    result = (matrix_[0][0] * matrix_[1][1]) - (matrix_[1][0] * matrix_[0][1]);
  } else {
    S21Matrix minor(rows_ - 1, cols_ - 1);
    for (int i = 0; i < cols_; i++) {
      NewMatrixForMinor(minor, 0, i);
      result += alg_add * matrix_[0][i] * minor.Determinant();
      alg_add = -alg_add;
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  ExceptionMatrixDoesntExist();
  ExceptionMatrixIsNotSquare();
  double determinant = Determinant();
  if (fabs(determinant) < EPS) {
    throw std::length_error("Matrix determinant is 0");
  }
  S21Matrix alg_matrix_transpose(rows_, cols_);
  alg_matrix_transpose = CalcComplements().Transpose();
  alg_matrix_transpose.MulNumber(1 / determinant);
  return alg_matrix_transpose;
}

void S21Matrix::OutOfRange(int i, int j) const {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::length_error("Out of range");
  }
}

void S21Matrix::ExceptionMatrixDoesntExist() const {
  if (matrix_ == nullptr) {
    throw std::length_error("Matrix doesn't exist");
  }
}

void S21Matrix::ExceptionMatrixDoesntExist(const S21Matrix& other) const {
  if (matrix_ == nullptr || other.matrix_ == nullptr) {
    throw std::length_error("Matrix doesn't exist");
  }
}

void S21Matrix::ExceptionMatrixIsNotSquare() const {
  if (rows_ != cols_) {
    throw std::length_error("Matrix is not square");
  }
}

void S21Matrix::ExceptionDifferentDimensionsForMul(
    const S21Matrix& other) const {
  if (cols_ != other.rows_) {
    throw std::length_error("Different dimensions for mul");
  }
}

void S21Matrix::ExceptionDifferentDimensions(const S21Matrix& other) const {
  if ((rows_ != other.rows_) || (cols_ != other.cols_)) {
    throw std::length_error("Different dimensions");
  }
}

void S21Matrix::NewMatrixForMinor(S21Matrix& minor, int crossed_out_rows,
                                  int crossed_out_column) noexcept {
  int minor_index_i = 0;
  int minor_index_j = 0;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (i != crossed_out_rows && j != crossed_out_column) {
        minor.matrix_[minor_index_i][minor_index_j++] = matrix_[i][j];
        if (minor_index_j == (cols_ - 1)) {
          minor_index_i++;
          minor_index_j = 0;
        }
      }
    }
  }
}

void S21Matrix::AllocNewMatrix(double**& matrix, int rows, int cols) noexcept {
  matrix = new double*[rows]();
  for (int i = 0; i < rows; i++) {
    matrix[i] = new double[cols]();
  }
  for (int i = 0; i < rows && i < rows_; i++) {
    for (int j = 0; j < cols; j++) {
      matrix[i][j] = matrix_[i][j];
    }
  }
}