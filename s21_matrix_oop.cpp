#include "s21_matrix_oop.h"

/**
 * @brief Default constructor for the S21Matrix class.
 * Initializes a matrix object with no elements (rows and columns set to 0) and
 * a null pointer for the matrix data.
 *
 * This constructor initializes an S21Matrix object without allocating any
 * memory for the matrix data, setting both dimensions (rows and columns) to 0
 * and the matrix pointer to nullptr. It's used when creating an empty matrix
 * object without specifying dimensions upfront.
 */
S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0), matrix_(nullptr) {}

/**
 * @brief Constructs an S21Matrix object with specified dimensions.
 *
 * This constructor creates an S21Matrix object with the specified number of
 * rows and columns. It dynamically allocates memory for the matrix data,
 * ensuring each element is initialized to zero.
 *
 * @param rows The number of rows for the matrix.
 * @param cols The number of columns for the matrix.
 *
 * @throws std::invalid_argument If either rows or cols is less than or equal to
 * zero, indicating invalid dimensions.
 */
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  assert(rows >= 0 && cols >= 0 && "cols/rows must be more then zero");
  this->matrix_ = new double*[rows];
  assert(this->matrix_ != nullptr && "memory not be alloceted");
  for (int i = 0; i < rows; i++) {
    this->matrix_[i] = new double[cols]{};
    assert(this->matrix_[i] != nullptr && "memory not be alloceted");
  }
}

/**
 * @brief Copy constructor for the S21Matrix class.
 *
 * This constructor creates an S21Matrix object based on the provided Matrix
 * object. It initializes the S21Matrix object with the dimensions of the Matrix
 * and fills it with the Matrix's contents.
 *
 * @param matrix The Matrix object to copy from.
 */
S21Matrix::S21Matrix(const Matrix& matrix)
    : S21Matrix(matrix.size(), matrix.begin()->size()) {
  is_square_matrix(matrix);
  FillingMatrix(matrix);
}

/**
 * @brief Constructor for the S21Matrix class that initializes the matrix from
 * an Array object.
 *
 * This constructor creates an S21Matrix object based on the provided Array
 * object, initializing it with the size of the array as the number of rows and
 * 1 as the number of columns. It then fills the newly constructed matrix with
 * the values from the array.
 *
 * @param array The Array object whose elements are used to fill the matrix. The
 * size of the array determines the number of rows in the matrix, while the
 * number of columns is set to 1.
 */
S21Matrix::S21Matrix(const Array& array) : S21Matrix(array.size(), 1) {
  FillingMatrix(array);
}

/**
 * @brief Destructor for the S21Matrix class.
 *
 * Cleans up the dynamically allocated memory associated with the matrix data.
 * If the matrix_ pointer is not null, it deletes each row of the matrix and
 * then deletes the array of pointers itself, setting matrix_ to nullptr and
 * resetting rows_ and cols_ to 0. Ensures that the S21Matrix object leaves no
 * dangling pointers or memory leaks upon destruction.
 */
S21Matrix::~S21Matrix() noexcept { Destructor(*this); }

/**
 * @brief Fills the S21Matrix object with the contents of another Matrix object.
 *
 * This method copies the elements from the provided Matrix object into the
 * current S21Matrix instance. It iterates over each element of the input Matrix
 * and assigns them to the corresponding positions in the S21Matrix's internal
 * matrix storage, effectively transferring the data from the input Matrix to
 * the S21Matrix.
 *
 * @param matrix A constant reference to the Matrix object whose contents are
 * copied into the S21Matrix. The dimensions of the input Matrix must match the
 * dimensions of the S21Matrix instance for correct operation.
 */
void S21Matrix::FillingMatrix(const Matrix& matrix) {
  int row = 0, col = 0;
  for (const auto& rows : matrix) {
    col = 0;
    for (const auto& elem : rows) {
      this->matrix_[row][col] = elem;
      col++;
    }
    row++;
  }
}

/**
 * @brief Fills the S21Matrix object with the contents of an Array object.
 *
 * This method initializes the S21Matrix object with the elements from the
 * provided Array, setting each element of the matrix to the corresponding value
 * from the array. The array is assumed to represent a single column of the
 * matrix, with its size determining the number of rows in the matrix, and all
 * elements are placed in the first column of the S21Matrix.
 *
 * @param array A constant reference to the Array object whose elements are used
 * to fill the matrix. The size of the array determines the number of rows in
 * the matrix, while the number of columns is implicitly set to 1, placing all
 * elements in the first column of the S21Matrix.
 */
void S21Matrix::FillingMatrix(const Array& array) {
  int i = 0;
  for (const auto& elem : array) {
    this->matrix_[i][0] = elem;
    i++;
  }
}

/**
 * @brief Updates the S21Matrix object with the contents of another Matrix
 * object.
 *
 * This method sets the contents of the current S21Matrix instance to match
 * those of the provided Matrix object. It first checks if the provided Matrix
 * is square and then asserts that the dimensions of the current matrix match
 * those of the input Matrix. After validating the dimensions, it fills the
 * current matrix with the contents of the input Matrix.
 *
 * @param matrix A constant reference to the Matrix object whose contents are
 * used to update the current S21Matrix instance. The dimensions of the input
 * Matrix must match the dimensions of the S21Matrix instance for correct
 * operation.
 */
void S21Matrix::set(const Matrix& matrix) {
  is_square_matrix(matrix);
  assert(this->rows_ == (int)matrix.size() && "rows incorrect");
  assert(this->cols_ == (int)matrix.begin()->size() && "rows incorrect");
  FillingMatrix(matrix);
}

/**
 * @brief Sets the S21Matrix object with the contents of an Array object.
 *
 * This method updates the S21Matrix object to match the contents of the
 * provided Array. It asserts that the number of rows in the S21Matrix matches
 * the size of the Array, ensuring the matrix can accommodate the new data.
 * Then, it fills the matrix with the elements from the Array, assuming the
 * Array represents a single column of data.
 *
 * @param array A constant reference to the Array object whose elements are used
 * to update the S21Matrix. The size of the Array determines the number of rows
 * in the matrix, and all elements are placed in the first column of the
 * S21Matrix.
 */
void S21Matrix::set(const Array& array) {
  assert(this->rows_ == int(array.size()) && "rows incorrect");
  FillingMatrix(array);
}

/**
 * @brief Returns a pointer to the matrix data.
 *
 * This method provides access to the internal representation of the matrix,
 * returning a pointer to the dynamically allocated 2D array that stores the
 * matrix elements. It allows direct manipulation of the matrix data but should
 * be used with caution to avoid unintended side effects.
 *
 * @return Pointer to the dynamically allocated 2D array representing the
 * matrix.
 */
double** S21Matrix::matrix() const noexcept { return this->matrix_; }

/**
 * @brief Returns the number of rows in the matrix.
 *
 * This method retrieves the number of rows in the S21Matrix object, providing a
 * way to query the dimensionality of the matrix along the row axis.
 *
 * @return The number of rows in the matrix.
 */
int S21Matrix::get_rows() const noexcept { return this->rows_; }

/**
 * @brief Returns the number of columns in the matrix.
 *
 * This method retrieves the number of columns in the S21Matrix object,
 * providing a way to query the dimensionality of the matrix along the column
 * axis.
 *
 * @return The number of columns in the matrix.
 */
int S21Matrix::get_cols() const noexcept { return this->cols_; }

/**
 * @brief Accesses the element at the specified position in the matrix.
 *
 * Provides access to the element located at the specified row and column in the
 * matrix. Checks are performed to ensure that the indices are valid, throwing
 * an exception if they are out of bounds.
 *
 * @param i The row index of the element to access.
 * @param j The column index of the element to access.
 * @return Reference to the element at the specified position in the matrix.
 * @throws std::out_of_range If the indices are out of bounds (less than zero or
 * greater than or equal to the number of rows/columns).
 */
double& S21Matrix::operator()(int i, int j) const {
  if (i < 0 || j < 0) {
    char msg[] = "Index of row or column could not be less zero";
    throw std::out_of_range(msg);
  }
  if (i >= this->rows_) {
    char msg[] = "Index of row must be less the number of matrix rows";
    throw std::out_of_range(msg);
  }
  if (j >= this->cols_) {
    char msg[] = "Index of column must be less the number of matrix columns";
    throw std::out_of_range(msg);
  }
  return matrix_[i][j];
}

/**
 * @brief Accesses the row at the specified index in the matrix.
 *
 * Provides access to the internal representation of the matrix, returning a
 * pointer to the dynamically allocated array representing the specified row. It
 * allows direct manipulation of the matrix data but should be used with caution
 * to avoid unintended side effects.
 *
 * @param index The index of the row to access.
 * @return Pointer to the dynamically allocated array representing the specified
 * row in the matrix.
 */
double* S21Matrix::operator[](int index) { return this->matrix_[index]; }

/**
 * @brief Checks if the given Matrix is square.
 *
 * This method verifies that the provided Matrix object is not null, has a
 * non-zero size, and that the number of rows equals the number of columns,
 * indicating a square matrix. It asserts that the matrix is not null, has valid
 * dimensions, and that each row has the same number of columns, ensuring the
 * matrix is rectangular.
 *
 * @param matrix A constant reference to the Matrix object to check. The method
 * checks if the matrix is square by comparing the number of rows to the number
 * of columns and ensuring consistency in the number of columns across all rows.
 */
void S21Matrix::is_square_matrix(const Matrix& matrix) {
  assert(matrix.size() > 0 && "Matrix is null!");

  const std::size_t rows = matrix.size();
  assert(rows != 0u && "Invalid rows in initializer matrix!");

  const std::size_t cols = matrix.begin()->size();
  assert(cols != 0u && "Invalid cols in initializer matrix!");
  const char msg[] =
      "Invalid cols in init-r matrix! Matrix is not rectangular!";
  for (const auto& line : matrix) assert((line.size() == cols) && msg);
}

/**
 * @brief Checks if two matrices are approximately equal within a small epsilon.
 *
 * This method compares two matrices for equality, considering a small epsilon
 * value to account for floating-point precision errors. It first checks if the
 * dimensions of both matrices match. If they do, it then iterates through each
 * element, comparing the absolute difference between corresponding elements of
 * the two matrices to a predefined epsilon value (1e-5). If all differences are
 * less than this epsilon, the matrices are considered equal.
 *
 * @param first The first matrix to compare.
 * @param second The second matrix to compare.
 * @return true if the matrices are approximately equal within the epsilon
 * value, false otherwise.
 */
bool S21Matrix::EqMatrix(const S21Matrix& first, const S21Matrix& second) {
  if (first.rows_ != second.rows_ || first.cols_ != second.cols_) return false;

  for (int i = 0; i < first.rows_; i++)
    for (int j = 0; j < first.cols_; j++)
      if (std::abs(first.matrix_[i][j] - second.matrix_[i][j]) > 1e-5)
        return false;
  return true;
}

/**
 * @brief Checks if this matrix is approximately equal to another S21Matrix
 * within a small epsilon.
 *
 * Compares this matrix with another S21Matrix for equality, considering a small
 * epsilon value to account for floating-point precision errors. It first checks
 * if the dimensions match and then iterates through each element, comparing the
 * absolute difference between corresponding elements to a predefined epsilon
 * value (1e-5). If all differences are less than this epsilon, the matrices are
 * considered equal.
 *
 * @param other The other S21Matrix to compare with.
 * @return true if the matrices are approximately equal within the epsilon
 * value, false otherwise.
 */
bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  return EqMatrix(*this, other);
}

/**
 * @brief Overloaded equality operator for comparing two S21Matrix objects.
 *
 * Compares this matrix with another S21Matrix for equality, considering a small
 * epsilon value to account for floating-point precision errors. It first checks
 * if the dimensions match and then iterates through each element, comparing the
 * absolute difference between corresponding elements to a predefined epsilon
 * value (1e-5). If all differences are less than this epsilon, the matrices are
 * considered equal.
 *
 * @param other The other S21Matrix to compare with.
 * @return true if the matrices are approximately equal within the epsilon
 * value, false otherwise.
 */
bool S21Matrix::operator==(const S21Matrix& other) const {
  return EqMatrix(*this, other);
}

/**
 * @brief Copy constructor for the S21Matrix class.
 *
 * This constructor creates an S21Matrix object based on the provided S21Matrix
 * object. It initializes the new S21Matrix object with the dimensions of the
 * other matrix and fills it with the contents of the other matrix.
 *
 * @param other The S21Matrix object to copy from.
 */
S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

/**
 * @brief Assignment operator for the S21Matrix class.
 *
 * Assigns the contents of another S21Matrix object to the current instance. If
 * the current instance is the same as the other, it returns immediately to
 * avoid self-assignment. Otherwise, it first deallocates the existing matrix
 * data, then copies the dimensions and contents from the other matrix to the
 * current instance.
 *
 * @param other The S21Matrix object to copy from.
 * @return Reference to the current instance after assignment.
 */
S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this == &other) return *this;
  Destructor(*this);

  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = new double*[this->rows_];
  for (int i = 0; i < this->rows_; i++) {
    this->matrix_[i] = new double[this->cols_];
    std::copy(other.matrix_[i], other.matrix_[i] + this->cols_,
              this->matrix_[i]);
  }
  return *this;
}

/**
 * @brief Assignment operator for assigning a Matrix object to the current
 * S21Matrix instance.
 *
 * Assigns the contents of a Matrix object to the current S21Matrix instance. It
 * creates a temporary S21Matrix object from the given Matrix, then moves the
 * contents of this temporary object into the current instance using move
 * semantics to efficiently transfer ownership of resources.
 *
 * @param matrix The Matrix object to assign from.
 * @return Reference to the current instance after assignment.
 */
S21Matrix S21Matrix::operator=(const Matrix& matrix) {
  S21Matrix copy_object(matrix);
  *this = std::move(copy_object);
  return *this;
}

/**
 * @brief Move constructor for the S21Matrix class.
 *
 * This constructor creates an S21Matrix object by moving from the provided
 * S21Matrix object. It initializes the new object with the dimensions and
 * matrix data of the other object, leaving the other object in a valid but
 * unspecified state by setting its matrix pointer to nullptr and its dimensions
 * to zero.
 *
 * @param other The S21Matrix object to move from.
 */
S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

/**
 * @brief Adds the contents of another S21Matrix object to the current instance.
 *
 * This method updates the current S21Matrix instance by adding the
 * corresponding elements of the provided S21Matrix object to it. It first
 * checks if the dimensions of both matrices match, throwing an exception if
 * they do not. If the dimensions match, it iterates through each element,
 * adding the values from the other matrix to the current matrix.
 *
 * @param other The S21Matrix object whose contents are added to the current
 * instance. The dimensions of the input S21Matrix must match the dimensions of
 * the current instance for correct operation.
 * @throws std::invalid_argument If the sizes of the matrices do not match.
 */
void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw std::invalid_argument("Sizes of matrices not match.");
  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++)
      this->matrix_[i][j] += other.matrix_[i][j];
}

/**
 * @brief Overloaded addition operator for adding two S21Matrix objects.
 *
 * Creates a new S21Matrix object that is the result of adding the current
 * instance to another S21Matrix object. It first creates a copy of the current
 * instance and then adds the other matrix to this copy using the SumMatrix
 * method. The resulting matrix is returned, leaving the original matrices
 * unchanged.
 *
 * @param other The S21Matrix object to be added to the current instance.
 * @return A new S21Matrix object that is the sum of the current instance and
 * the other matrix.
 */
S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix tmp = *this;
  tmp.SumMatrix(other);
  return tmp;
}

/**
 * @brief Overloaded addition assignment operator for adding another S21Matrix
 * object to the current instance.
 *
 * Adds the contents of another S21Matrix object to the current instance in
 * place, modifying the current instance to be the sum of the original matrix
 * and the other matrix. It uses the SumMatrix method to perform the addition.
 *
 * @param other The S21Matrix object to be added to the current instance.
 * @return Reference to the current instance after addition.
 */
S21Matrix S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return *this;
}

/**
 * @brief Subtracts the contents of another S21Matrix object from the current
 * instance.
 *
 * This method updates the current S21Matrix instance by subtracting the
 * corresponding elements of the provided S21Matrix object from it. It first
 * checks if the dimensions of both matrices match, throwing an exception if
 * they do not. If the dimensions match, it iterates through each element,
 * subtracting the values from the other matrix from the current matrix.
 *
 * @param other The S21Matrix object whose contents are subtracted from the
 * current instance. The dimensions of the input S21Matrix must match the
 * dimensions of the current instance for correct operation.
 * @throws std::invalid_argument If the sizes of the matrices do not match.
 */
void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (this->rows_ != other.rows_ || this->cols_ != other.cols_)
    throw std::invalid_argument("Sizes of matrices not match.");
  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++)
      this->matrix_[i][j] -= other.matrix_[i][j];
}

/**
 * @brief Overloaded subtraction operator for subtracting two S21Matrix objects.
 *
 * Creates a new S21Matrix object that is the result of subtracting another
 * S21Matrix object from the current instance. It first creates a copy of the
 * current instance and then subtracts the other matrix from this copy using the
 * SubMatrix method. The resulting matrix is returned, leaving the original
 * matrices unchanged.
 *
 * @param other The S21Matrix object to be subtracted from the current instance.
 * @return A new S21Matrix object that is the difference of the current instance
 * and the other matrix.
 */
S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix tmp = *this;
  tmp.SubMatrix(other);
  return tmp;
}

/**
 * @brief Overloaded subtraction assignment operator for subtracting another
 * S21Matrix object from the current instance.
 *
 * Subtracts the contents of another S21Matrix object from the current instance
 * in place, modifying the current instance to be the difference of the original
 * matrix and the other matrix. It uses the SubMatrix method to perform the
 * subtraction.
 *
 * @param other The S21Matrix object to be subtracted from the current instance.
 * @return Reference to the current instance after subtraction.
 */
S21Matrix S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return *this;
}

/**
 * @brief Multiplies the current S21Matrix instance by a scalar.
 *
 * This method updates the current S21Matrix instance by multiplying each of its
 * elements by a given number. It iterates through each element of the matrix,
 * applying the multiplication in place.
 *
 * @param num The scalar by which to multiply the matrix elements.
 */
void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++) this->matrix_[i][j] *= num;
}

/**
 * @brief Multiplies the current S21Matrix instance by another S21Matrix object.
 *
 * This method updates the current S21Matrix instance to be the product of its
 * multiplication with another S21Matrix object. It first checks if the number
 * of columns of the current matrix matches the number of rows of the other
 * matrix, throwing an exception if they do not. It then performs matrix
 * multiplication, creating a new result matrix with dimensions matching the
 * current matrix's rows and the other matrix's columns. The result is assigned
 * to the current instance using move semantics.
 *
 * @param other The S21Matrix object to multiply with the current instance. The
 * number of columns of the current matrix must equal the number of rows of the
 * other matrix for the multiplication to be valid.
 * @throws std::invalid_argument If the number of columns of the current matrix
 * does not match the number of rows of the other matrix.
 */
void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (this->cols_ != other.rows_) {
    throw std::invalid_argument(
        "Number of columns of the matrix has to be equal number or rows of "
        "second matrix");
  }
  int rows = this->rows_;
  int columns = other.cols_;
  S21Matrix result(this->rows_, other.cols_);

  for (int i = 0; i != rows; ++i) {
    for (int j = 0; j != columns; ++j) {
      for (int k = 0; k != this->cols_; ++k) {
        result(i, j) += (*this)(i, k) * other(k, j);
      }
    }
  }
  *this = std::move(result);
}

/**
 * @brief Multiplies the current S21Matrix instance by another S21Matrix object
 * in place.
 *
 * This operator modifies the current instance to be the product of its
 * multiplication with another S21Matrix object. It internally calls the
 * MulMatrix function to perform the matrix multiplication and then returns a
 * reference to the current instance.
 *
 * @param other The S21Matrix object to multiply with the current instance. The
 * number of columns of the current matrix must equal the number of rows of the
 * other matrix for the multiplication to be valid.
 * @return Reference to the current instance after multiplication.
 */
S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

/**
 * @brief Multiplies the current S21Matrix instance by a scalar value in place.
 *
 * This operator modifies the current instance by multiplying each of its
 * elements by a given scalar value. It internally calls the MulNumber function
 * to perform the multiplication and then returns a reference to the current
 * instance.
 *
 * @param value The scalar value by which to multiply the matrix elements.
 * @return Reference to the current instance after multiplication.
 */
S21Matrix& S21Matrix::operator*=(const double& value) {
  MulNumber(value);
  return *this;
}

/**
 * @brief Creates a new S21Matrix object that is the result of multiplying the
 * current instance by another S21Matrix object.
 *
 * This operator does not modify the current instance but instead creates a
 * temporary copy, performs the multiplication using the MulMatrix function, and
 * returns the resulting matrix.
 *
 * @param other The S21Matrix object to multiply with the current instance. The
 * number of columns of the current matrix must equal the number of rows of the
 * other matrix for the multiplication to be valid.
 * @return A new S21Matrix object that is the product of the current instance
 * and the other matrix.
 */
S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix tmp = *this;
  tmp.MulMatrix(other);
  return tmp;
}

/**
 * @brief Creates a new S21Matrix object that is the result of multiplying the
 * current instance by a scalar value.
 *
 * This operator does not modify the current instance but instead creates a
 * temporary copy, applies the scalar multiplication using the MulNumber
 * function, and returns the resulting matrix.
 *
 * @param value The scalar value by which to multiply the matrix elements.
 * @return A new S21Matrix object that is the product of the current instance
 * and the scalar value.
 */
S21Matrix S21Matrix::operator*(const double& value) const {
  S21Matrix res{*this};
  res.MulNumber(value);
  return res;
}
/**
 * @brief Creates a new S21Matrix object that is the result of multiplying a
 * scalar value by another S21Matrix object.
 *
 * This operator allows for multiplication of a scalar value with a matrix from
 * the left, effectively calling the right-hand side operator* implementation.
 *
 * @param value The scalar value to multiply with the matrix elements.
 * @param matrix The S21Matrix object to be multiplied by the scalar value.
 * @return A new S21Matrix object that is the product of the scalar value and
 * the matrix.
 */
S21Matrix operator*(const double& value, const S21Matrix& matrix) {
  return matrix * value;
}

/**
 * @brief Creates a new S21Matrix object that is the transpose of the current
 * instance.
 *
 * This method generates a new matrix where the rows and columns of the original
 * matrix are swapped. Each element at position (n, m) in the original matrix is
 * placed at position (m, n) in the new matrix.
 *
 * @return A new S21Matrix object that is the transpose of the current instance.
 */
S21Matrix S21Matrix::Transpose() const {
  S21Matrix new_matrix(this->cols_, this->rows_);
  for (int n = 0; n < this->rows_; n++) {
    for (int m = 0; m < this->cols_; m++) {
      new_matrix[m][n] = this->matrix_[n][m];
    }
  }
  return new_matrix;
}

/**
 * @brief Calculates the determinant of the current S21Matrix instance.
 *
 * The determinant is a special number that can be calculated from a square
 * matrix. This method computes the determinant of the matrix. The determinant
 * helps us find the inverse of a matrix, tells us things about the matrix that
 * are useful in systems of linear equations, calculus and more.
 *
 * @throws std::logic_error if the matrix is not square (number of rows is not
 * equal to the number of columns) or if the matrix is empty (number of rows or
 * columns is zero).
 * @return The determinant of the matrix as a double value.
 */
double S21Matrix::Determinant() const {
  char msg[] =
      "Determinant could be computed only for not empty square matrices";
  if (this->get_rows() != this->get_cols() || this->get_rows() == 0)
    throw std::logic_error(msg);

  int size = this->get_rows();
  if (size == 1) return (*this)(0, 0);

  double determinant = 0;
  for (int j = 0, sign = 1; j != size; ++j) {
    S21Matrix minor_matrix = this->MinorMatrix(0, j);
    double minor_determinant = minor_matrix.Determinant();
    determinant += sign * (*this)(0, j) * minor_determinant;
    sign *= -1;
  }
  return determinant;
}

/**
 * @brief Calculates the minor matrix for a given row and column of the current
 * S21Matrix instance.
 *
 * The minor matrix is a smaller matrix derived from the original matrix by
 * removing a specific row and column. It is used in the calculation of
 * determinants and inverses of matrices.
 *
 * @param row The row to be removed from the original matrix.
 * @param column The column to be removed from the original matrix.
 * @throws std::logic_error if the matrix is not square or if it is empty or has
 * only one row/column.
 * @throws std::range_error if the specified row or column is out of range.
 * @return A new S21Matrix object representing the minor matrix.
 */
S21Matrix S21Matrix::MinorMatrix(int row, int column) const {
  if (this->get_rows() != this->get_cols()) {
    throw std::logic_error(
        "Minor matrix could be computed only for square matrices");
  }
  if (this->get_rows() <= 1 || this->get_cols() <= 1) {
    throw std::logic_error(
        "The empty or one sized matrices doesn't have minors");
  }
  if (row < 0 || row >= this->get_rows()) {
    throw std::range_error("The row number out of range for the matrix.");
  }
  if (column < 0 || column >= this->get_cols()) {
    throw std::range_error("The columns number out of range for the matrix.");
  }

  int size = this->get_rows() - 1;
  S21Matrix minor(size, size);

  for (int i = 0, minor_i = 0; i != this->get_rows(); ++i) {
    if (i == row) continue;
    for (int j = 0, minor_j = 0; j < this->get_cols(); ++j) {
      if (j == column) continue;
      minor(minor_i, minor_j) = (*this)(i, j);
      ++minor_j;
    }
    ++minor_i;
  }
  return minor;
}

/**
 * @brief Calculates the matrix of cofactors for the current S21Matrix instance.
 *
 * The matrix of cofactors is used in the calculation of the inverse of a
 * matrix. Each element of the matrix is replaced by the determinant of its
 * minor multiplied by (-1)^(i+j), where i and j are the row and column indices
 * of the element.
 *
 * @throws std::logic_error if the matrix is not square or empty.
 * @return A new S21Matrix object representing the matrix of cofactors.
 */
S21Matrix S21Matrix::CalcComplements() const {
  if (this->get_rows() != this->get_cols() || this->get_rows() == 0) {
    throw std::logic_error(
        "Complements matrix could be computed only for not empty square "
        "matrices");
  }
  int size = this->get_rows();
  S21Matrix complements(size, size);

  for (int i = 0; i != size; ++i) {
    for (int j = 0; j != size; ++j) {
      int sign = ((i + j) % 2 == 0) ? 1 : -1;
      complements(i, j) = sign * (*this).MinorMatrix(i, j).Determinant();
    }
  }
  return complements;
}

/**
 * @brief Calculates the inverse of the current S21Matrix instance.
 *
 * The inverse of a matrix is a matrix that, when multiplied by the original
 * matrix, yields the identity matrix. This method computes the inverse using
 * the matrix of cofactors and the determinant.
 *
 * @throws std::logic_error if the matrix is not square or empty.
 * @throws std::overflow_error if the determinant of the matrix is zero,
 * indicating that the matrix is singular (non-invertible).
 * @return A new S21Matrix object representing the inverse of the current
 * instance, or an empty matrix if the inverse does not exist.
 */
S21Matrix S21Matrix::InverseMatrix() const {
  if (this->get_rows() != this->get_cols() || this->get_rows() == 0) {
    throw std::logic_error(
        "Inverse matrix could be computed only for not empty square "
        "matrices");
  }
  double determinant = this->Determinant();
  if (fabs(determinant) < 1e-6) {
    throw std::overflow_error("The determinant is equal zero");
  }
  int size = this->get_rows();
  S21Matrix transposed_complements = (*this).CalcComplements().Transpose();
  S21Matrix inverse_matrix(size, size);
  for (int i = 0; i != size; ++i) {
    for (int j = 0; j != size; ++j) {
      inverse_matrix(i, j) = transposed_complements(i, j) / determinant;
    }
  }
  return inverse_matrix;
}

/**
 * @brief Cleans up the dynamically allocated memory associated with the matrix
 * data of an S21Matrix object.
 *
 * This function is designed to be used as a helper for the destructor and other
 * member functions that need to deallocate the matrix data of an S21Matrix
 * object safely. It checks if the matrix_ pointer is not null, then deletes
 * each row of the matrix and the array of pointers itself, setting matrix_ to
 * nullptr and resetting rows_ and cols_ to 0. This ensures that the S21Matrix
 * object leaves no dangling pointers or memory leaks upon destruction.
 *
 * @param matrix The S21Matrix object whose dynamically allocated memory is to
 * be cleaned up.
 */
void Destructor(S21Matrix& matrix) noexcept {
  if (matrix.matrix_ == nullptr) return;
  for (int i = 0; i < matrix.rows_; i++) delete[] matrix.matrix_[i];
  delete[] matrix.matrix_;
  matrix.matrix_ = nullptr;
  matrix.cols_ = 0;
  matrix.rows_ = 0;
}

/**
 * @brief Changes the size of the current S21Matrix instance to the specified
 * dimensions.
 *
 * This method adjusts the size of the current S21Matrix instance to the
 * specified number of rows and columns. If the new size is different from the
 * current size, it creates a new matrix of the specified size, copies the
 * existing elements to the new matrix (if they fit), and then replaces the old
 * matrix with the new one. If the new size is the same as the current size, the
 * method does nothing.
 *
 * @param rows The new number of rows for the matrix.
 * @param cols The new number of columns for the matrix.
 */
void S21Matrix::ChangeSize(int rows, int cols) noexcept {
  if (rows == this->rows_ && cols == this->cols_) return;
  S21Matrix new_matrix(rows, cols);
  int x = (rows < this->rows_) ? rows : this->rows_;
  int y = (cols < this->cols_) ? cols : this->cols_;
  for (int i = 0; i < x; i++)
    for (int j = 0; j < y; j++) new_matrix[i][j] = this->matrix_[i][j];
  Destructor(*this);
  *this = new_matrix;
}

/**
 * @brief Deletes a specified row from the current S21Matrix instance.
 *
 * This method removes the specified row from the matrix, adjusting the
 * matrix's dimensions accordingly. It creates a new matrix with one less row,
 * copies the remaining elements from the original matrix, and then replaces
 * the original matrix with this new matrix. The operation is not performed if
 * the specified row index is out of bounds.
 *
 * @param row_index The index of the row to be deleted from the matrix.
 * @throws std::invalid_argument if the row index is less than 0 or greater
 * than or equal to the current number of rows in the matrix.
 */
void S21Matrix::del_row(int row_index) {
  if (row_index < 0 || row_index > this->rows_)
    throw std::invalid_argument("row_index incorrect");

  S21Matrix new_matrix(this->rows_ - 1, this->cols_);
  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++) {
      // skip deleting row_index
      if (i == row_index) continue;
      new_matrix[i][j] = this->matrix_[i][j];
    }
  Destructor(*this);
  *this = new_matrix;
}

/**
 * @brief Deletes a specified column from the current S21Matrix instance.
 *
 * This method removes the specified column from the matrix, adjusting the
 * matrix's dimensions accordingly. It creates a new matrix with one less
 * column, copies the remaining elements from the original matrix, and then
 * replaces the original matrix with this new matrix. The operation is not
 * performed if the specified column index is out of bounds.
 *
 * @param col_index The index of the column to be deleted from the matrix.
 * @throws std::invalid_argument if the column index is less than 0 or greater
 * than or equal to the current number of columns in the matrix.
 */
void S21Matrix::del_col(int col_index) {
  if (col_index < 0 || col_index > this->cols_)
    throw std::invalid_argument("col_index incorrect");

  S21Matrix new_matrix(this->rows_, this->cols_ - 1);
  for (int i = 0; i < this->rows_; i++)
    for (int j = 0; j < this->cols_; j++) {
      // skip deleting col_index
      if (j == col_index) continue;
      new_matrix[i][j] = this->matrix_[i][j];
    }
  Destructor(*this);
  *this = new_matrix;
}

/**
 * @brief Inserts a new row into the current S21Matrix instance at the
 * specified index.
 *
 * This method inserts a new row from the provided Array into the matrix at
 * the specified row index. If the row index is not provided or is -1, the new
 * row is inserted at the end of the matrix. It creates a new matrix with one
 * more row, copies the existing elements and the new row into the new matrix,
 * and then replaces the original matrix with this new matrix. The operation
 * is not performed if the specified row index is out of bounds.
 *
 * @param arr The Array object containing the elements to insert as a new row.
 * The size of the array must match the number of columns in the matrix.
 * @param row_index The index at which the new row should be inserted. If -1,
 * the row is inserted at the end. If the index is out of bounds, an exception
 * is thrown.
 * @throws std::invalid_argument if the row index is less than 0 or greater
 * than the current number of columns in the matrix, or if the size of the
 * array does not match the number of columns.
 */
void S21Matrix::insert_row(Array arr, int row_index) {
  if (row_index < 0 || row_index > this->cols_)
    throw std::invalid_argument("row_index incorrect");
  S21Matrix new_matrix(this->rows_ + 1, this->cols_);

  double new_arr[arr.size()];
  int i = 0;
  for (auto&& elem : arr) {
    new_arr[i] = elem;
    i++;
  }

  for (int i = 0; i < this->rows_ + 1; i++)
    for (int j = 0; j < this->cols_; j++)
      if (i == row_index) {
        new_matrix[i][j] = new_arr[j];
      } else {
        new_matrix[i][j] = this->matrix_[i][j];
      }
  Destructor(*this);
  *this = new_matrix;
}

/**
 * @brief Inserts a new column into the current S21Matrix instance at the
 * specified index.
 *
 * This method inserts a new column from the provided Array into the matrix at
 * the specified column index. If the column index is not provided or is -1,
 * the new column is inserted at the end of the matrix. It creates a new
 * matrix with one more column, copies the existing elements and the new
 * column into the new matrix, and then replaces the original matrix with this
 * new matrix. The operation is not performed if the specified column index is
 * out of bounds.
 *
 * @param arr The Array object containing the elements to insert as a new
 * column. The size of the array must match the number of rows in the matrix.
 * @param col_index The index at which the new column should be inserted. If
 * -1, the column is inserted at the end. If the index is out of bounds, an
 * exception is thrown.
 * @throws std::invalid_argument if the column index is less than 0 or greater
 * than the current number of columns in the matrix, or if the size of the
 * array does not match the number of rows.
 */
void S21Matrix::insert_col(Array arr, int col_index) {
  if (col_index < 0 || col_index > this->cols_)
    throw std::invalid_argument("col_index incorrect");

  if (int(arr.size()) > this->rows_) throw std::invalid_argument("invalid col");
  S21Matrix new_matrix(this->rows_, this->cols_ + 1);

  double new_arr[arr.size()];
  int i = 0;
  for (auto&& elem : arr) {
    new_arr[i] = elem;
    i++;
  }

  for (int i = 0; i < this->rows_; i++) {
    for (int j = 0; j < this->cols_ + 1; j++) {
      if (j == col_index) {
        new_matrix.matrix_[i][j] = new_arr[i];  // insert new row
      } else {
        if (j < col_index) {  // copy elem from old matrix
          new_matrix.matrix_[i][j] = this->matrix_[i][j];
        } else {
          new_matrix.matrix_[i][j] = this->matrix_[i][j - 1];
        }
      }
    }
  }

  Destructor(*this);
  *this = new_matrix;
}

/**
 * @brief Adds a row to the end of the matrix.
 *
 * This method adds a new row, represented by an Array object, to the end of the
 * matrix. The size of the new row must match the number of columns in the
 * matrix. If the sizes do not match, an exception will be thrown.
 *
 * @param arr An Array object containing the elements of the new row to be
 * added.
 * @throws std::invalid_argument If the size of arr does not match the number of
 * columns in the matrix.
 */
void S21Matrix::add_row(Array arr) { return insert_row(arr, this->rows_); }

/**
 * @brief Adds a column to the end of the matrix.
 *
 * This method adds a new column, represented by an Array object, to the end of
 * the matrix. The size of the new column must match the number of rows in the
 * matrix. If the sizes do not match, an exception will be thrown.
 *
 * @param arr An Array object containing the elements of the new column to be
 * added.
 * @throws std::invalid_argument If the size of arr does not match the number of
 * rows in the matrix.
 */
void S21Matrix::add_col(Array arr) { return insert_col(arr, this->cols_); }