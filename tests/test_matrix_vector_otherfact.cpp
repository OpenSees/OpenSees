#include "DummyStream.h"
#include "Matrix.h"
#include "Vector.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

DummyStream testError;
OPS_Stream *opserrPtr = &testError;

namespace {

bool nearlyEqual(double lhs, double rhs)
{
  const double scale = 1.0 + std::max(std::fabs(lhs), std::fabs(rhs));
  return std::fabs(lhs - rhs) <= 1.0e-12 * scale;
}

void require(bool condition, const char *message)
{
  if (!condition) {
    std::cerr << "FAIL: " << message << "\n";
    std::exit(EXIT_FAILURE);
  }
}

Matrix makeMatrix(int rows, int cols, double offset)
{
  Matrix result(rows, cols);
  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
      result(i, j) = offset + 1.0 + 2.0 * i + 3.0 * j;
  return result;
}

Vector makeVector(int size, double offset)
{
  Vector result(size);
  for (int i = 0; i < size; ++i)
    result(i) = offset + 1.0 + 2.0 * i;
  return result;
}

Matrix scaled(const Matrix &value, double factor)
{
  Matrix result(value);
  for (int j = 0; j < result.noCols(); ++j)
    for (int i = 0; i < result.noRows(); ++i)
      result(i, j) *= factor;
  return result;
}

Matrix add(const Matrix &lhs, const Matrix &rhs, double rhsFactor)
{
  Matrix result(lhs);
  for (int j = 0; j < result.noCols(); ++j)
    for (int i = 0; i < result.noRows(); ++i)
      result(i, j) += rhs(i, j) * rhsFactor;
  return result;
}

Matrix multiply(const Matrix &lhs, const Matrix &rhs)
{
  Matrix result(lhs.noRows(), rhs.noCols());
  for (int j = 0; j < rhs.noCols(); ++j)
    for (int i = 0; i < lhs.noRows(); ++i)
      for (int k = 0; k < lhs.noCols(); ++k)
        result(i, j) += lhs(i, k) * rhs(k, j);
  return result;
}

Matrix transpose(const Matrix &value)
{
  Matrix result(value.noCols(), value.noRows());
  for (int j = 0; j < value.noCols(); ++j)
    for (int i = 0; i < value.noRows(); ++i)
      result(j, i) = value(i, j);
  return result;
}

void requireMatrix(const Matrix &actual, const Matrix &expected, const char *message)
{
  require(actual.noRows() == expected.noRows() && actual.noCols() == expected.noCols(), message);
  for (int j = 0; j < actual.noCols(); ++j)
    for (int i = 0; i < actual.noRows(); ++i)
      require(nearlyEqual(actual(i, j), expected(i, j)), message);
}

void requireVector(const Vector &actual, const Vector &expected, const char *message)
{
  require(actual.Size() == expected.Size(), message);
  for (int i = 0; i < actual.Size(); ++i)
    require(nearlyEqual(actual(i), expected(i)), message);
}

} // namespace

int main()
{
  const Matrix A = makeMatrix(2, 2, 0.0);
  const Matrix B = makeMatrix(2, 2, 10.0);

  Matrix actual = A;
  require(actual.addMatrix(1.0, B, 0.0) == 0, "addMatrix zero factor return");
  requireMatrix(actual, A, "addMatrix zero factor preserves this");
  actual = A;
  actual.addMatrix(2.0, B, 0.0);
  requireMatrix(actual, scaled(A, 2.0), "addMatrix zero factor scales this");
  actual = A;
  actual.addMatrix(1.0, B, -1.0);
  requireMatrix(actual, add(A, B, -1.0), "addMatrix negative-one factor");
  actual = A;
  actual.addMatrix(2.0, B, 2.5);
  requireMatrix(actual, add(scaled(A, 2.0), B, 2.5), "addMatrix general factors");

  actual = A;
  actual.addMatrixTranspose(1.0, B, 0.0);
  requireMatrix(actual, A, "addMatrixTranspose zero factor");
  actual = A;
  actual.addMatrixTranspose(1.0, B, -1.0);
  requireMatrix(actual, add(A, transpose(B), -1.0), "addMatrixTranspose negative-one factor");
  actual = A;
  actual.addMatrixTranspose(2.0, B, 2.5);
  requireMatrix(actual, add(scaled(A, 2.0), transpose(B), 2.5), "addMatrixTranspose general factors");

  const Matrix B23 = makeMatrix(2, 3, 1.0);
  const Matrix C32 = makeMatrix(3, 2, 2.0);
  const Matrix product = multiply(B23, C32);
  actual = A;
  actual.addMatrixProduct(1.0, B23, C32, 0.0);
  requireMatrix(actual, A, "addMatrixProduct zero factor");
  actual = A;
  actual.addMatrixProduct(0.0, B23, C32, 1.0);
  requireMatrix(actual, product, "addMatrixProduct positive-one factor");
  actual = A;
  actual.addMatrixProduct(1.0, B23, C32, -1.0);
  requireMatrix(actual, add(A, product, -1.0), "addMatrixProduct negative-one factor");
  actual = A;
  actual.addMatrixProduct(2.0, B23, C32, 2.5);
  requireMatrix(actual, add(scaled(A, 2.0), product, 2.5), "addMatrixProduct general factors");

  const Matrix A32 = makeMatrix(3, 2, 3.0);
  const Matrix B32 = makeMatrix(3, 2, 4.0);
  const Matrix transposeProduct = multiply(transpose(A32), B32);
  actual = A;
  actual.addMatrixTransposeProduct(1.0, A32, B32, 0.0);
  requireMatrix(actual, A, "addMatrixTransposeProduct zero factor");
  actual = A;
  actual.addMatrixTransposeProduct(1.0, A32, B32, -1.0);
  requireMatrix(actual, add(A, transposeProduct, -1.0), "addMatrixTransposeProduct negative-one factor");
  actual = A;
  actual.addMatrixTransposeProduct(2.0, A32, B32, 2.5);
  requireMatrix(actual, add(scaled(A, 2.0), transposeProduct, 2.5), "addMatrixTransposeProduct general factors");

  const Matrix T32 = makeMatrix(3, 2, 5.0);
  const Matrix B33 = makeMatrix(3, 3, 6.0);
  const Matrix triple = multiply(transpose(T32), multiply(B33, T32));
  actual = A;
  actual.addMatrixTripleProduct(1.0, T32, B33, 0.0);
  requireMatrix(actual, A, "addMatrixTripleProduct zero factor");
  actual = A;
  actual.addMatrixTripleProduct(1.0, T32, B33, -1.0);
  requireMatrix(actual, add(A, triple, -1.0), "addMatrixTripleProduct negative-one factor");
  actual = A;
  actual.addMatrixTripleProduct(2.0, T32, B33, 2.5);
  requireMatrix(actual, add(scaled(A, 2.0), triple, 2.5), "addMatrixTripleProduct general factors");

  const Matrix C32b = makeMatrix(3, 2, 7.0);
  const Matrix tripleGeneral = multiply(transpose(A32), multiply(B33, C32b));
  actual = A;
  actual.addMatrixTripleProduct(1.0, A32, B33, C32b, 0.0);
  requireMatrix(actual, A, "general triple product zero factor");
  actual = A;
  actual.addMatrixTripleProduct(1.0, A32, B33, C32b, -1.0);
  requireMatrix(actual, add(A, tripleGeneral, -1.0), "general triple product negative-one factor");
  actual = A;
  actual.addMatrixTripleProduct(2.0, A32, B33, C32b, 2.5);
  requireMatrix(actual, add(scaled(A, 2.0), tripleGeneral, 2.5), "general triple product general factors");

  const Vector v = makeVector(2, 8.0);
  const Vector w = makeVector(2, 9.0);
  Vector vector = v;
  vector.addVector(1.0, w, 0.0);
  requireVector(vector, v, "addVector zero factor");
  vector = v;
  vector.addVector(0.0, w, -1.0);
  Vector expectedVector = w;
  expectedVector *= -1.0;
  requireVector(vector, expectedVector, "addVector negative-one factor");
  vector = v;
  vector.addVector(2.0, w, 2.5);
  for (int i = 0; i < expectedVector.Size(); ++i)
    expectedVector(i) = 2.0 * v(i) + 2.5 * w(i);
  requireVector(vector, expectedVector, "addVector general factors");

  const Matrix M23 = makeMatrix(2, 3, 11.0);
  const Vector x3 = makeVector(3, 12.0);
  Vector matrixVector = v;
  matrixVector.addMatrixVector(1.0, M23, x3, 0.0);
  requireVector(matrixVector, v, "addMatrixVector zero factor");
  Vector Mx = makeVector(2, 0.0);
  for (int i = 0; i < M23.noRows(); ++i) {
    Mx(i) = 0.0;
    for (int j = 0; j < M23.noCols(); ++j)
      Mx(i) += M23(i, j) * x3(j);
  }
  matrixVector = v;
  matrixVector.addMatrixVector(0.0, M23, x3, -1.0);
  Mx *= -1.0;
  requireVector(matrixVector, Mx, "addMatrixVector negative-one factor");
  Mx *= -1.0;
  matrixVector = v;
  matrixVector.addMatrixVector(2.0, M23, x3, 2.5);
  for (int i = 0; i < expectedVector.Size(); ++i)
    expectedVector(i) = 2.0 * v(i) + 2.5 * Mx(i);
  requireVector(matrixVector, expectedVector, "addMatrixVector general factors");

  const Matrix M32 = makeMatrix(3, 2, 13.0);
  const Vector x3b = makeVector(3, 14.0);
  Vector transposeVector = v;
  transposeVector.addMatrixTransposeVector(1.0, M32, x3b, 0.0);
  requireVector(transposeVector, v, "addMatrixTransposeVector zero factor");
  Vector Mtx = makeVector(2, 0.0);
  for (int j = 0; j < M32.noCols(); ++j) {
    Mtx(j) = 0.0;
    for (int i = 0; i < M32.noRows(); ++i)
      Mtx(j) += M32(i, j) * x3b(i);
  }
  transposeVector = v;
  transposeVector.addMatrixTransposeVector(0.0, M32, x3b, -1.0);
  Mtx *= -1.0;
  requireVector(transposeVector, Mtx, "addMatrixTransposeVector negative-one factor");
  Mtx *= -1.0;
  transposeVector = v;
  transposeVector.addMatrixTransposeVector(2.0, M32, x3b, 2.5);
  for (int i = 0; i < expectedVector.Size(); ++i)
    expectedVector(i) = 2.0 * v(i) + 2.5 * Mtx(i);
  requireVector(transposeVector, expectedVector, "addMatrixTransposeVector general factors");

  std::cout << "matrix/vector otherFact regression tests passed\n";
  return EXIT_SUCCESS;
}
