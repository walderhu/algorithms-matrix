#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iostream>

using Array = std::initializer_list<double>;
using Matrix = std::initializer_list<Array>;

/**
 * @class S21Matrix
 * @brief A class representing a matrix with dynamic size and operations for
matrix manipulations.
 *
 * This class provides functionalities for creating and manipulating matrices of
dynamic size, including basic operations such as addition, subtraction,
 * multiplication, and finding the determinant, inverse, and transpose of a
matrix. It supports operations with both matrices and initializer lists for
 * initialization and comparison. The class ensures that matrices are always in
a valid state, throwing exceptions for invalid operations like non-square
 * matrix determinants or inverses.
 * matrix determinants or inverses.
 *
 * The S21Matrix class is designed to work with matrices of double values,
offering a flexible interface for matrix arithmetic and transformations. It
 * includes constructors for creating matrices from scratch, copying from other
matrices, or initializing from initializer lists. It also provides a
 * destructor to clean up allocated memory.
 *
 * Operations include setting matrix values, comparing matrices, arithmetic
operations (addition, subtraction, multiplication), and advanced operations
 * like transposition, determinant calculation, minor matrix calculation,
cofactor matrix calculation, and inverse matrix calculation. It also supports
 * scalar multiplication and division.
 *
 * The class uses dynamic memory allocation for storing matrix elements,
ensuring that matrices can be of any size, limited only by available memory.
 * It includes checks for matrix validity (e.g., square matrices for
determinants and inverses) and throws exceptions for invalid operations.
 *
 * Example usage:
 * @code
 * S21Matrix m1(3, 3); // Creates a 3x3 matrix
 * m1.set({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}); // Sets values using an initializer
list
 * S21Matrix m2 = m1.Transpose(); // Creates a transpose of m1
 * double det = m1.Determinant(); // Calculates the determinant of m1
 * @endcode
 *
 * The class is designed to be efficient and flexible, allowing for both direct
manipulation of matrix elements and high-level operations. It uses modern
 */
class S21Matrix final {
 private:
  int rows_, cols_;
  double** matrix_;

 public:
  S21Matrix() noexcept;
  S21Matrix(int rows, int cols);
  S21Matrix(const Matrix& matrix);
  S21Matrix(const Array& array);
  ~S21Matrix() noexcept;

  void FillingMatrix(const Matrix& matrix);
  void FillingMatrix(const Array& array);
  void print() const;
  void set(const Matrix& matrix);
  void set(const Array& array);
  double** matrix() const noexcept;
  int get_rows() const noexcept;
  int get_cols() const noexcept;

  double& operator()(int i, int j) const;
  double* operator[](int index);
  static void is_square_matrix(const Matrix& matrix);

  static bool EqMatrix(const S21Matrix& first, const S21Matrix& second);
  bool EqMatrix(const S21Matrix& other) const;
  bool operator==(const S21Matrix& other) const;
  bool operator==(const Matrix matrix) const;
  bool operator==(const Array array) const;
  bool EqMatrix(double** arr) const;

  S21Matrix(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix operator=(const Matrix& matrix);
  S21Matrix operator=(const Array& array);
  S21Matrix(S21Matrix&& other);

  void SumMatrix(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator+=(const S21Matrix& other);

  void SubMatrix(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator-=(const S21Matrix& other);

  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double& value);
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(const double& value) const;
  friend S21Matrix operator*(const double& value, const S21Matrix& matrix);

  S21Matrix Transpose() const;
  double Determinant() const;
  S21Matrix MinorMatrix(int row, int column) const;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;

  friend void Destructor(S21Matrix& matrix) noexcept;
  void ChangeSize(int rows, int cols) noexcept;
  void del_row(int row_index);
  void del_col(int col_index);

  void insert_row(Array arr, int row_index);
  void insert_col(Array arr, int col_index);
  void add_row(Array arr);
  void add_col(Array arr);
};

#endif
