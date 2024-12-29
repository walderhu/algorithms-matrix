#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(AddTests, NullMatrix) {
  S21Matrix m;
  int rows = m.get_rows();
  int cols = m.get_cols();
  double **matrix = m.matrix();
  bool cond = rows == 0 && cols == 0 && matrix == nullptr;
  EXPECT_TRUE(cond);
}

TEST(AddTests, EmptyMatrixStock) {
  S21Matrix m(3, 3);
  for (int i = 0; i < m.get_rows(); i++)
    for (int j = 0; j < m.get_cols(); j++) EXPECT_TRUE(m[i][j] == 0);
}

TEST(AddTests, EmptyMatrixFill) {
  S21Matrix m({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
  for (int i = 0; i < m.get_rows(); i++)
    for (int j = 0; j < m.get_cols(); j++) EXPECT_TRUE(m[i][j] == 0);
}

TEST(AddTests, IsSquareMatrixCorrect) {
  Matrix m = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  try {
    S21Matrix::is_square_matrix(m);
    EXPECT_TRUE(true);
  } catch (const std::exception &e) {
    EXPECT_TRUE(false);
  }
}

TEST(AddTests, staples) {
  Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m(matrix);
  int arr[3][3] = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) EXPECT_TRUE(m[i][j] == arr[i][j]);
}

TEST(AddTests, Set) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  Matrix matrix2 = {
      {4, 0, 3},
      {1, 0, 9},
      {7, 8, -6},
  };
  int arr1[3][3] = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  int arr2[3][3] = {
      {4, 0, 3},
      {1, 0, 9},
      {7, 8, -6},
  };

  S21Matrix m(matrix1);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) EXPECT_TRUE(m[i][j] == arr1[i][j]);

  m.set(matrix2);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) EXPECT_TRUE(m[i][j] == arr2[i][j]);
}

TEST(AddTests, StaticEqMatrixTrue) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  Matrix matrix2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(matrix2);
  bool cond = S21Matrix::EqMatrix(m1, m2);
  EXPECT_TRUE(cond);
}

TEST(AddTests, StaticEqMatrixFalse) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
  };
  Matrix matrix2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(matrix2);
  bool cond = S21Matrix::EqMatrix(m1, m2);
  EXPECT_FALSE(cond);
}

TEST(AddTests, EqMatrixFalse) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
  };
  Matrix matrix2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(matrix2);
  bool cond = m1.EqMatrix(m2);
  EXPECT_FALSE(cond);
}

TEST(AddTests, EqMatrixTrue) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  Matrix matrix2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(matrix2);
  bool cond = m1.EqMatrix(m2);
  EXPECT_TRUE(cond);
}

TEST(AddTests, operatorEqMatrixFalse) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
  };
  Matrix matrix2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(matrix2);
  bool cond = m1 == m2;
  EXPECT_FALSE(cond);
}

TEST(AddTests, operatorEqMatrixTrue) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  Matrix matrix2 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(matrix2);
  bool cond = m1 == m2;
  EXPECT_TRUE(cond);
}
TEST(AddTests, copy_constructor) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2(m1);
  bool cond = m1 == m2;
  EXPECT_TRUE(cond);
}

TEST(AddTests, operator_transfer_constructor) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2;
  m2 = m1;
  bool cond = m1 == m2;
  EXPECT_TRUE(cond);
}

TEST(AddTests, operator_transfer_true) {
  Matrix matrix1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  S21Matrix m1(matrix1);
  S21Matrix m2 = matrix1;
  bool cond = m1 == m2;
  EXPECT_TRUE(cond);
}

TEST(AddTests, SomeTest1) {
  S21Matrix m1 = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  S21Matrix m2 = {
      {4, 6, 3},
      {4, 1, 6},
      {7, 4, 9},
  };
  bool cond = m1 == m2;
  EXPECT_FALSE(cond);
}

TEST(MatrixEq, EmptyMatricesShouldBeEqual) {
  S21Matrix first;
  S21Matrix second;
  bool result = first.EqMatrix(second);
  EXPECT_TRUE(result);
}

TEST(MatrixEq, NotEqualIfRowsDiffers) {
  S21Matrix first(1, 1);
  S21Matrix second(2, 1);
  bool result = first.EqMatrix(second);
  EXPECT_FALSE(result);
}

TEST(MatrixEq, NotEqualIfColumnsDiffers) {
  S21Matrix first(1, 1);
  S21Matrix second(1, 2);
  bool result = first.EqMatrix(second);
  EXPECT_FALSE(result);
}

TEST(MatrixEq, EqualsIfValuesDiffersButNotThanEpsilon) {
  S21Matrix first(1, 1);
  first[0][0] = 2.00000001;
  S21Matrix second(1, 1);
  second[0][0] = 2.00000003;
  bool result = first.EqMatrix(second);
  EXPECT_TRUE(result);
}

TEST(SumMatrix, SummingEmptyMatrixReturnEmptyMatrix) {
  S21Matrix first_empty;
  const S21Matrix second_empty;
  first_empty.SumMatrix(second_empty);
  EXPECT_TRUE(first_empty == second_empty);
}

TEST(SumMatrix, SumMatrixActuallySumMatrices) {
  S21Matrix first = {1, 2, 3, 4, 5, 6};
  const S21Matrix second = {1, 2, 3, 4, 5, 6};
  first.SumMatrix(second);
  EXPECT_TRUE(first.EqMatrix({2, 4, 6, 8, 10, 12}));
}

TEST(SumMatrix, ThrowIfRowsNotMatch) {
  S21Matrix first = {1, 2, 3, 4, 5, 6};
  const S21Matrix second = {1, 2};
  EXPECT_THROW(first.SumMatrix(second), std::invalid_argument);
}

TEST(SumMatrix, ThrowIfColumnsNotMatch) {
  S21Matrix first = {1, 2, 3, 4, 5, 6};
  const S21Matrix second = {1, 2, 3};
  EXPECT_THROW(first.SumMatrix(second), std::invalid_argument);
}

TEST(SumMatrix, SumMatrixOperatorWorksOk) {
  const S21Matrix first = {1, 2};
  const S21Matrix second = {1, 2};
  const S21Matrix expected = {2, 4};
  S21Matrix result = first + second;
  EXPECT_TRUE(result == expected);
}

TEST(SumMatrix, SumSelfMatrixOperatorWorksOk) {
  S21Matrix matrix = {1, 2};
  const S21Matrix second = {1, 2};
  const S21Matrix expected = {2, 4};
  matrix += second;
  EXPECT_TRUE(matrix == expected);
}

TEST(SubtractMatrix, SubtractingEmptyMatricesDoNothing) {
  S21Matrix first_empty;
  const S21Matrix second_empty;
  first_empty.SubMatrix(second_empty);
  EXPECT_TRUE(first_empty == second_empty);
}

TEST(SubtractMatrix, ActuallySubtractingMatrices) {
  S21Matrix first = {
      {1, 2, 3},
      {4, 5, 6},
  };
  S21Matrix second = {
      {1, 2, 3},
      {4, 5, 6},
  };
  first.SubMatrix(second);
  S21Matrix result(2, 3);
  bool cond = first == result;
  EXPECT_TRUE(cond);
}

TEST(SubtractMatrix, ThrowIfRowsNotMatch) {
  S21Matrix first = {1, 2, 3, 4, 5, 6};
  const S21Matrix second = {1, 2};

  EXPECT_THROW(first.SubMatrix(second), std::invalid_argument);
}

TEST(SubtractMatrix, ThrowIfColumnsNotMatch) {
  S21Matrix first = {1, 2, 3, 4, 5, 6};
  const S21Matrix second = {1, 2, 3};
  EXPECT_THROW(first.SubMatrix(second), std::invalid_argument);
}

TEST(SubtractMatrix, SubtractingMatrixOperatorWorksOk) {
  const S21Matrix first = {1, 2};
  const S21Matrix second = {1, 2};
  S21Matrix result = first - second;

  S21Matrix correct = {0, 0};
  EXPECT_TRUE(result == correct);
}

TEST(SubtractMatrix, SubtractingSelfMatrixOperatorWorksOk) {
  S21Matrix matrix = {1, 2};
  const S21Matrix second = {1, 2};
  const S21Matrix expected(1, 2);
  matrix -= second;
  EXPECT_FALSE(matrix == expected);
}

TEST(MultiplyByNumber, MultiplyingEmptyMatrixDoNothing) {
  S21Matrix matrix;
  matrix.MulNumber(10);
  EXPECT_EQ(matrix.get_rows(), 0);
  EXPECT_EQ(matrix.get_cols(), 0);
  EXPECT_EQ(matrix, S21Matrix());
}

TEST(MultiplyByNumber, MultiplyingByNumberMultiplyAllCells) {
  S21Matrix matrix = {
      {1, 2, 3, 4},
      {5, 6, 7, 8},
  };

  matrix.MulNumber(2);

  EXPECT_EQ(matrix.get_rows(), 2);
  EXPECT_EQ(matrix.get_cols(), 4);
  S21Matrix result = {
      {2, 4, 6, 8},
      {10, 12, 14, 16},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(MultiplyByNumber, MultiplyingByNumberOperator) {
  const S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };
  S21Matrix result = matrix * 2;

  EXPECT_EQ(result.get_rows(), 2);
  EXPECT_EQ(result.get_cols(), 3);
  EXPECT_TRUE(result == S21Matrix({
                            {2, 4, 6},
                            {8, 10, 12},
                        }));
}

TEST(MultiplyByNumber, MultiplyingByNumberOperatorAlsoWorks) {
  const S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };
  S21Matrix result = 2 * matrix;

  EXPECT_EQ(result.get_rows(), 2);
  EXPECT_EQ(result.get_cols(), 3);
  EXPECT_TRUE(result == S21Matrix({
                            {2, 4, 6},
                            {8, 10, 12},
                        }));
}

TEST(MultiplyByNumber, SelfMultiplyingByNumberOperatorWorks) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };
  matrix *= 2;
  EXPECT_EQ(matrix.get_rows(), 2);
  EXPECT_EQ(matrix.get_cols(), 3);
  EXPECT_TRUE(matrix == S21Matrix({
                            {2, 4, 6},
                            {8, 10, 12},
                        }));
}

TEST(MultiplyByNumber, EmptyMatrixCouldBeMultiplied) {
  S21Matrix empty_matrix;
  S21Matrix result = empty_matrix * 10;
  EXPECT_TRUE(result == empty_matrix);
}

TEST(MulMatrix, EmptyCompatibleMatricesReturnEmptyMatrix) {
  S21Matrix first(3, 2);
  S21Matrix second(2, 3);
  first.MulMatrix(second);
  EXPECT_EQ(first.get_rows(), 3);
  EXPECT_EQ(first.get_cols(), 3);
  EXPECT_TRUE(first == S21Matrix({
                           {0, 0, 0},
                           {0, 0, 0},
                           {0, 0, 0},
                       }));
}

TEST(MulMatrix, MulAllRowsAndColumns) {
  S21Matrix first = {
      {1, 2},
      {3, 4},
      {5, 6},
  };
  S21Matrix second = {
      {1, 2, 3},
      {4, 5, 6},
  };
  first.MulMatrix(second);
  EXPECT_EQ(first.get_rows(), 3);
  EXPECT_EQ(first.get_cols(), 3);
  EXPECT_TRUE(first == S21Matrix({
                           {9, 12, 15},
                           {19, 26, 33},
                           {29, 40, 51},
                       }));
}

TEST(MulMatrix, ThrowIfMatricesNotCompatibleForMul) {
  S21Matrix first(1, 2);
  S21Matrix second(3, 2);
  EXPECT_THROW(first.MulMatrix(second), std::invalid_argument);
}

TEST(MulMatrix, MulMatrixOperator) {
  S21Matrix first(2, 2);
  first = {
      {1, 2},
      {3, 4},
  };
  S21Matrix second(2, 1);
  Array arr1 = {
      {1},
      {2},
  };
  second.set(arr1);

  S21Matrix expected(2, 1);
  Array arr2 = {
      {5},
      {11},
  };
  expected.set(arr2);
  S21Matrix result = first * second;
  EXPECT_TRUE(result == expected);
}

TEST(MulMatrix, MulSelfMatrixOperator) {
  S21Matrix first(2, 2);
  first = {
      {1, 2},
      {3, 4},
  };
  S21Matrix second(2, 1);
  Array arr1 = {
      {1},
      {2},
  };
  second.set(arr1);

  S21Matrix expected(2, 1);
  Array arr2 = {
      {5},
      {11},
  };
  expected.set(arr2);
  first *= second;
  EXPECT_TRUE(first == expected);
}

TEST(MatrixTranspose, EmptyMatrixTransposedToItself) {
  const S21Matrix empty_matrix;
  S21Matrix transposed = empty_matrix.Transpose();
  EXPECT_TRUE(transposed == empty_matrix);
}

TEST(MatrixTranspose, MatrixActuallyTransposed) {
  //           1 4   1 2 3
  // A = A^T = 2 5 = 4 5 6
  //           3 6
  const S21Matrix matrix = {
      {1, 4},
      {2, 5},
      {3, 6},
  };
  S21Matrix transposed = matrix.Transpose();
  EXPECT_EQ(transposed.get_rows(), 2);
  EXPECT_EQ(transposed.get_cols(), 3);
  EXPECT_TRUE(transposed == S21Matrix({
                                {1, 2, 3},
                                {4, 5, 6},
                            }));
}

TEST(MatrixDeterminant, MatrixDeterminantComputedOk) {
  const S21Matrix matrix = {
      {1, 2, 3, 4},
      {5, 6, 7, 8},
      {9, 10, 11, 12},
      {13, 14, 15, 16},
  };

  double determinant = matrix.Determinant();

  EXPECT_NEAR(determinant, 0.0, 1e-6);
}

TEST(MatrixDeterminant, MatrixDeterminantForOneElementMatrix) {
  const S21Matrix matrix = {686.686595};

  double determinant = matrix.Determinant();

  EXPECT_NEAR(determinant, 686.686595, 1e-6);
}

TEST(MatrixDeterminant, MatrixSizeSeverComputesOkAsWell) {
  const S21Matrix matrix = {
      {3.54155, 2.53027, 2.52268, 3.32609, 1.74077, 1.84826, 2.31548},
      {1.79850, 1.03137, 3.21930, 2.23851, 3.69014, 1.86757, 2.58197},
      {3.12373, 2.50464, 2.81140, 3.95159, 3.93592, 2.55369, 2.85939},
      {2.53200, 2.17887, 3.89360, 2.57050, 1.35048, 2.86216, 2.85716},
      {2.45656, 2.44057, 3.60225, 2.02151, 2.66313, 3.54608, 1.99800},
      {2.73763, 3.62892, 3.16649, 1.46655, 1.63051, 3.30205, 1.16198},
      {2.99739, 1.10405, 3.75781, 1.69789, 2.66463, 2.54331, 1.13451},
  };

  double determinant = matrix.Determinant();

  EXPECT_NEAR(determinant, 24.111673073047285806, 1e-6);
}

TEST(MatrixDeterminant, ThrowIfMatrixNotSquare) {
  const S21Matrix matrix(2, 3);

  EXPECT_THROW(matrix.Determinant(), std::logic_error);
}

TEST(MatrixDeterminant, ThrowIfMatrixIsEmpty) {
  const S21Matrix empty_matrix;
  EXPECT_THROW(empty_matrix.Determinant(), std::logic_error);
}

class MinorMatrixTest : public ::testing::Test {
 protected:
  const S21Matrix matrix4x4;

  MinorMatrixTest()
      : matrix4x4({
            {1, 2, 3, 4},
            {5, 6, 7, 8},
            {9, 10, 11, 12},
            {13, 14, 15, 16},
        }) {}
};

TEST_F(MinorMatrixTest, MinorMatrixByFirstRowAndFirstColumnCreated) {
  S21Matrix minor = matrix4x4.MinorMatrix(0, 0);

  EXPECT_EQ(minor.get_rows(), 3);
  EXPECT_EQ(minor.get_cols(), 3);
  EXPECT_TRUE(minor == S21Matrix({
                           {6, 7, 8},
                           {10, 11, 12},
                           {14, 15, 16},
                       }));
}

TEST_F(MinorMatrixTest, MinorMatrixByMiddleRowAndColumnCreated) {
  S21Matrix minor = matrix4x4.MinorMatrix(2, 1);

  EXPECT_EQ(minor.get_rows(), 3);
  EXPECT_EQ(minor.get_cols(), 3);
  EXPECT_TRUE(minor == S21Matrix({
                           {1, 3, 4},
                           {5, 7, 8},
                           {13, 15, 16},
                       }));
}

TEST_F(MinorMatrixTest, ThrowIfMatrixNotSquare) {
  S21Matrix matrix(3, 4);

  EXPECT_THROW(matrix.MinorMatrix(0, 0), std::logic_error);
}

TEST_F(MinorMatrixTest, ThrowIfMatrixEmpty) {
  S21Matrix matrix;

  EXPECT_THROW(matrix.MinorMatrix(0, 0), std::logic_error);
}

TEST_F(MinorMatrixTest, ThrowIfMatrixSizeOne) {
  S21Matrix matrix(1, 1);

  EXPECT_THROW(matrix.MinorMatrix(0, 0), std::logic_error);
}

TEST_F(MinorMatrixTest, ThrowIfRowsOutOfRange) {
  EXPECT_THROW(matrix4x4.MinorMatrix(4, 0), std::range_error);
}

TEST_F(MinorMatrixTest, ThrowIfColumnsOutOfRange) {
  EXPECT_THROW(matrix4x4.MinorMatrix(0, 4), std::range_error);
}

TEST(MatrixCalcComplements, ComplementFor3x3MatrixCalculated) {
  const S21Matrix matrix = {
      {1., 2., 3.},
      {3., 2., 1.},
      {7., 5., 2.},
  };

  S21Matrix complements = matrix.CalcComplements();

  EXPECT_EQ(complements.get_rows(), 3);
  EXPECT_EQ(complements.get_cols(), 3);

  S21Matrix m = {
      {-1., 1., 1.},
      {11., -19., 9.},
      {-4., 8., -4.},
  };
  EXPECT_TRUE(complements == m);
}

TEST(MatrixCalcComplements, ComplementsForNotSimpleMatrixCalculatedOk) {
  const S21Matrix matrix = {
      {1.80377, 3.93870, 3.13429, 2.28155},
      {1.39307, 1.05586, 2.21357, 2.20440},
      {2.74323, 2.41325, 3.86805, 2.73013},
      {2.29065, 3.09765, 1.84139, 3.86339},
  };

  S21Matrix complements = matrix.CalcComplements();

  S21Matrix m = {
      {-8.0642664633, 4.1987149757, 3.1661056480, -0.0941589509},
      {-19.1443430067, -4.4198224214, 8.7731760020, 10.7131854857},
      {15.1040957594, -1.3457695400, -1.9412358558, -6.9511236616},
      {5.0123523428, 0.9933255993, -5.5038169258, 1.8833757880},
  };
  EXPECT_TRUE(complements == m);
}

TEST(MatrixCalcComplements, ThrowIfMatrixIsEmpty) {
  S21Matrix empty_matrix;

  EXPECT_THROW(empty_matrix.CalcComplements(), std::logic_error);
}

TEST(MatrixCalcComplements, ThrowIfMatrixNotSquare) {
  S21Matrix matrix(1, 2);

  EXPECT_THROW(matrix.CalcComplements(), std::logic_error);
}

TEST(MatrixInverse, Regular3x3Case) {
  const S21Matrix matrix = {
      {1., 2., 3.},
      {3., 2., 1.},
      {7., 5., 2.},
  };
  const S21Matrix expected_result = {
      {-0.25, 2.75, -1.},
      {0.25, -4.75, 2.},
      {0.25, 2.25, -1.},
  };

  S21Matrix inverse_matrix = matrix.InverseMatrix();

  EXPECT_TRUE(inverse_matrix == expected_result);
}

TEST(MatrixInverse, TestResultForMatrix4x4) {
  const S21Matrix matrix = {
      {1.80377, 3.93870, 3.13429, 2.28155},
      {1.39307, 1.05586, 2.21357, 2.20440},
      {2.74323, 2.41325, 3.86805, 2.73013},
      {2.29065, 3.09765, 1.84139, 3.86339},
  };
  const S21Matrix expected_result = {
      {-0.6892499107, -1.6362600080, 1.2909415507, 0.4284039249},
      {0.3588626362, -0.3777606089, -0.1150224313, 0.0848991764},
      {0.2706058939, 0.7498401502, -0.1659167199, -0.4704092234},
      {-0.0080477312, 0.9156520525, -0.5941100018, 0.1609714411}};

  S21Matrix inverse_matrix = matrix.InverseMatrix();

  EXPECT_TRUE(inverse_matrix == expected_result);
}

TEST(MatrixInverse, ThrowIfDeterminantEqualsZero) {
  const S21Matrix matrix = {
      {1.80377, 3.93870, 3.13429, 2.28155},
      {1.39307, 1.05586, 2.21357, 2.20440},
      {2.74323, 2.41325, 3.86805, 2.73013},
      {2.74323, 2.41325, 3.86805, 2.73013},
  };

  EXPECT_THROW(matrix.InverseMatrix(), std::overflow_error);
}

TEST(MatrixInverse, ThrowIfMatrixEmpty) {
  const S21Matrix empty_matrix;

  EXPECT_THROW(empty_matrix.InverseMatrix(), std::logic_error);
}

TEST(MatrixInverse, ThrowIfMatrixNotSquare) {
  const S21Matrix matrix(1, 2);

  EXPECT_THROW(matrix.InverseMatrix(), std::logic_error);
}

TEST(MatrixOperators, GetItemValueByIndexes) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 2), 3);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 2), 6);
}

TEST(MatrixOperators, GetItemValueByIndexesForConst) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  EXPECT_EQ(matrix(0, 0), 1);
  EXPECT_EQ(matrix(0, 2), 3);
  EXPECT_EQ(matrix(1, 0), 4);
  EXPECT_EQ(matrix(1, 2), 6);
}

TEST(MatrixOperators, ThrowGetByIndexIfIndexOutOfRange) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  EXPECT_THROW(matrix(-1, 0), std::out_of_range);
  EXPECT_THROW(matrix(3, 0), std::out_of_range);
  EXPECT_THROW(matrix(0, -1), std::out_of_range);
  EXPECT_THROW(matrix(0, 3), std::out_of_range);
}

TEST(MatrixOperators, ThrowGetByIndexIfEmptyMatrix) {
  const S21Matrix matrix;

  EXPECT_THROW(matrix(0, 0), std::out_of_range);
}

TEST(AddTests, ChangeMatrix1) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  matrix.ChangeSize(2, 2);

  S21Matrix result = {
      {1, 2},
      {4, 5},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, ChangeMatrix2) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  matrix.ChangeSize(3, 3);

  S21Matrix result = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, ChangeMatrix3) {
  S21Matrix matrix = {
      {1, 2},
      {4, 5},
  };
  matrix.ChangeSize(3, 3);

  S21Matrix result = {
      {1, 2, 0},
      {4, 5, 0},
      {0, 0, 0},
  };

  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, RemoveRowCol) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  matrix.del_col(2);
  matrix.del_row(2);

  S21Matrix result = {
      {1, 2},
      {4, 5},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, RemoveRow) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  matrix.del_row(2);

  S21Matrix result = {
      {1, 2, 3},
      {4, 5, 6},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, RemoveCol) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  matrix.del_col(2);

  S21Matrix result = {
      {1, 2},
      {4, 5},
      {7, 8},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, AddRowCol) {
  S21Matrix matrix = {
      {1, 2},
      {4, 5},
  };
  matrix.add_row({7, 8});
  matrix.add_col({3, 6, 9});

  S21Matrix result = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, add_col) {
  S21Matrix matrix = {
      {1, 2},
      {4, 5},
      {7, 8},
  };
  matrix.add_col({3, 6, 9});
  S21Matrix result = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  EXPECT_TRUE(matrix == result);
}
TEST(AddTests, insert_col) {
  S21Matrix matrix = {
      {1, 2},
      {4, 5},
      {7, 8},
  };
  matrix.insert_col({3, 6, 9}, 1);
  S21Matrix result = {
      {1, 3, 2},
      {4, 6, 5},
      {7, 9, 8},
  };
  EXPECT_TRUE(matrix == result);
}

TEST(AddTests, ad_row) {
  S21Matrix matrix = {
      {1, 2, 3},
      {4, 5, 6},
  };

  matrix.add_row({7, 8, 9});

  S21Matrix result = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };

  EXPECT_TRUE(matrix == result);
}

TEST(S21MatrixTest, MoveConstructor) {
  // Создаем исходный объект матрицы
  S21Matrix matrix1 = {
      {1, 1, 1},
      {2, 2, 2},
      {3, 3, 3},
  };  // Перемещаем объект
  S21Matrix matrix2(std::move(matrix1));

  // Проверяем, что matrix2 правильно инициализирована
  EXPECT_EQ(matrix2(0, 0), 1.0);
  EXPECT_EQ(matrix2(1, 1), 2.0);
  EXPECT_EQ(matrix2(2, 2), 3.0);
  EXPECT_EQ(matrix2.get_rows(), 3);  // Проверяем количество строк
  EXPECT_EQ(matrix2.get_cols(), 3);  // Проверяем количество столбцов

  // Проверяем, что matrix1 теперь "пустая"
  EXPECT_EQ(matrix1.get_rows(), 0);  // После перемещения, rows_ должен быть 0
  EXPECT_EQ(matrix1.get_cols(), 0);  // После перемещения, cols_ должен быть 0
  EXPECT_EQ(matrix1.matrix(), nullptr);
}

TEST(AddTests, throwCheck) {
  S21Matrix matrix1 = {
      {1, 1, 1},
      {2, 2, 2},
      {3, 3, 3},
  };  // Перемещаем объект
  EXPECT_THROW(matrix1.del_col(55), std::invalid_argument);
  EXPECT_THROW(matrix1.del_row(55), std::invalid_argument);
  EXPECT_THROW(matrix1.insert_col({1, 2, 3}, 55), std::invalid_argument);
  EXPECT_THROW(matrix1.insert_row({1, 2, 3}, 55), std::invalid_argument);
}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}