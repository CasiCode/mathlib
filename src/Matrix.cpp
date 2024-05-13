/*
 * matrix.cpp
 */

#include "../include/Matrix.hpp"
#include <cmath>
#include <tuple>

#include <iostream> // DELETE LATER

#define EPS 1e-10

using std::domain_error;
using std::endl;
using std::istream;
using std::ostream;

/* PUBLIC MEMBER FUNCTIONS
 ********************************/

Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  allocSpace();
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] = 0;
    }
  }
}

Matrix::Matrix(double **a, int rows, int cols) : rows_(rows), cols_(cols) {
  allocSpace();
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] = a[i][j];
    }
  }
}

Matrix::Matrix() : rows_(1), cols_(1) {
  allocSpace();
  p[0][0] = 0;
}

Matrix::~Matrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] p[i];
  }
  delete[] p;
}

Matrix::Matrix(const Matrix &m) : rows_(m.rows_), cols_(m.cols_) {
  allocSpace();
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] = m.p[i][j];
    }
  }
}

Matrix &Matrix::operator=(const Matrix &m) {
  if (this == &m) {
    return *this;
  }

  if (rows_ != m.rows_ || cols_ != m.cols_) {
    for (int i = 0; i < rows_; ++i) {
      delete[] p[i];
    }
    delete[] p;

    rows_ = m.rows_;
    cols_ = m.cols_;
    allocSpace();
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] = m.p[i][j];
    }
  }
  return *this;
}

Matrix &Matrix::operator+=(const Matrix &m) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] += m.p[i][j];
    }
  }
  return *this;
}

Matrix &Matrix::operator-=(const Matrix &m) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] -= m.p[i][j];
    }
  }
  return *this;
}

Matrix &Matrix::operator*=(const Matrix &m) {
  Matrix temp(rows_, m.cols_);
  for (int i = 0; i < temp.rows_; ++i) {
    for (int j = 0; j < temp.cols_; ++j) {
      for (int k = 0; k < cols_; ++k) {
        temp.p[i][j] += (p[i][k] * m.p[k][j]);
      }
    }
  }
  return (*this = temp);
}

Matrix &Matrix::operator*=(double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] *= num;
    }
  }
  return *this;
}

Matrix &Matrix::operator/=(double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      p[i][j] /= num;
    }
  }
  return *this;
}

Matrix Matrix::operator^(int num) {
  Matrix temp(*this);
  return expHelper(temp, num);
}

void Matrix::swapRows(int r1, int r2) {
  double *temp = p[r1];
  p[r1] = p[r2];
  p[r2] = temp;
}

Matrix Matrix::transpose() {
  Matrix ret(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      ret.p[j][i] = p[i][j];
    }
  }
  return ret;
}

/* STATIC CLASS FUNCTIONS
 ********************************/

Matrix Matrix::createIdentity(int size) {
  Matrix temp(size, size);
  for (int i = 0; i < temp.rows_; ++i) {
    for (int j = 0; j < temp.cols_; ++j) {
      if (i == j) {
        temp.p[i][j] = 1;
      } else {
        temp.p[i][j] = 0;
      }
    }
  }
  return temp;
}

/*void QRdecomp(int n, int m, Matrix A, Matrix Q) { // TODO
  double r;
  double c;
  double s;
  double temp;
  int min;
  Q = Matrix::createIdentity(m);
  min = (m < n ? m : n);
  for (int i = 0; i < min; i++)
    for (int j = i + 1; j < m; j++)
      if (A(j, i) != 0) {
        r = hypot(A(i, i), A(j, i));
        c = static_cast<double>(A(i, i)) / r;
        s = static_cast<double>(A(j, i)) / r;
        for (int k = 0; k < n; k++) {
          temp = A(i, k);
          A(i, k) = c * A(i, k) + s * A(j, k);
          A(j, k) = -s * temp + c * A(j, k);
        }
        for (int k = 0; k < m; k++) {
          temp = Q(k, i);
          Q(k, i) = c * Q(k, i) + s * Q(k, j);
          Q(k, j) = -s * temp + c * Q(k, j);
        }
      }
}*/

/*Matrix Matrix::solve(Matrix A, Matrix b) {
  // Gaussian elimination
  for (int i = 0; i < A.rows_; ++i) {
    if (A.p[i][i] == 0) {
      // pivot 0 - throw error
      throw domain_error("Error: the coefficient matrix has 0 as a pivot. "
                         "Please fix the input and try again.");
    }
    for (int j = i + 1; j < A.rows_; ++j) {
      for (int k = i + 1; k < A.cols_; ++k) {
        A.p[j][k] -= A.p[i][k] * (A.p[j][i] / A.p[i][i]);
        if (A.p[j][k] < EPS && A.p[j][k] > -1 * EPS)
          A.p[j][k] = 0;
      }
      b.p[j][0] -= b.p[i][0] * (A.p[j][i] / A.p[i][i]);
      if (A.p[j][0] < EPS && A.p[j][0] > -1 * EPS)
        A.p[j][0] = 0;
      A.p[j][i] = 0;
    }
  }

  // Back substitution
  Matrix x(b.rows_, 1);
  x.p[x.rows_ - 1][0] = b.p[x.rows_ - 1][0] / A.p[x.rows_ - 1][x.rows_ - 1];
  if (x.p[x.rows_ - 1][0] < EPS && x.p[x.rows_ - 1][0] > -1 * EPS)
    x.p[x.rows_ - 1][0] = 0;
  for (int i = x.rows_ - 2; i >= 0; --i) {
    int sum = 0;
    for (int j = i + 1; j < x.rows_; ++j) {
      sum += A.p[i][j] * x.p[j][0];
    }
    x.p[i][0] = (b.p[i][0] - sum) / A.p[i][i];
    if (x.p[i][0] < EPS && x.p[i][0] > -1 * EPS)
      x.p[i][0] = 0;
  }

  return x;
}*/

Matrix Matrix::solve(Matrix A, Matrix b) {
  // Gaussian elimination
  for (int i = 0; i < A.cols_; ++i) { // rows_ ???
    int pivotRow = i;
    for (int j = i + 1; j < A.rows_; ++j) {
      if (std::fabs(A(i, j)) > std::fabs(A(pivotRow, i))) {
        pivotRow = j;
      }
    }
    if (pivotRow != i) {
      A.swapRows(i, pivotRow);
      b.swapRows(i, pivotRow);
    }
    if (A(i, i) == 0) {
      // pivot 0 - throw error
      throw domain_error("Error: the coefficient matrix has 0 as a pivot. "
                         "Please fix the input and try again.");
    }

    for (int j = i + 1; j < A.rows_; ++j) {
      for (int k = i + 1; k < A.cols_; ++k) {
        A(j, k) -= A(i, k) * (A(j, i) / A(i, i));
        if (A(j, k) < EPS && A(j, k) > -1 * EPS)
          A(j, k) = 0;
      }
      b(j, 0) -= b(i, 0) * (A(j, i) / A(i, i));
      if (b(j, 0) < EPS &&
          b(j, 0) > -1 * EPS) // changed A.p to b.p everywhere!!!
        b(j, 0) = 0;          // here too
      A(j, i) = 0;
    }
  }

  // Back substitution
  Matrix x(b.rows_, 1);

  x(x.rows_ - 1, 0) = b(x.rows_ - 1, 0) / A(x.rows_ - 1, x.rows_ - 1);
  if (x(x.rows_ - 1, 0) < EPS && x(x.rows_ - 1, 0) > -1 * EPS)
    x(x.rows_ - 1, 0) = 0;
  for (int i = x.rows_ - 2; i >= 0; --i) {
    int sum = 0;
    for (int j = i + 1; j < x.rows_; ++j) {
      sum += A(i, j) * x(j, 0);
    }
    x(i, 0) = (b(i, 0) - sum) / A(i, i);
    if (x(i, 0) < EPS && x(i, 0) > -1 * EPS)
      x(i, 0) = 0;
  }

  return x;
}

std::tuple<Matrix, Matrix, int, int>
Matrix::jacobiEigenvalue(int n, Matrix a, int maxIterations) {
  //  Purpose:
  //    JACOBI_EIGENVALUE computes the eigenvalues and eigenvectors of a
  //    real symmetric matrix, using Rutishauser's modfications of the classical
  //    Jacobi rotation method with threshold pivoting.
  //  Input:
  //    int N, the order of the matrix.
  //    Matrix A[N*N], the matrix, which must be square, real,
  //    and symmetric.
  //    int IT_MAX, the maximum number of iterations.
  //  Output:
  //    double V[N*N], the matrix of eigenvectors.
  //    double D[N], the eigenvalues, in descending order.
  //    int &IT_NUM, the total number of iterations.
  //    int &ROT_NUM, the total number of rotations.
  double *bw;
  double c;
  double g;
  double gapq;
  double h;
  double s;
  double t;
  double tau;
  double term;
  double termp;
  double termq;
  double theta;
  double thresh;
  double w;
  double *zw;

  Matrix eigenvecs = Matrix::createIdentity(n);
  Matrix eigenvals = Matrix(1, n);

  int iterations = 0;
  int rotations = 0;

  eigenvals = a.getMainDiagonal();

  bw = new double[n];
  zw = new double[n];

  for (int i = 0; i < n; i++) {
    bw[i] = eigenvals[i];
    zw[i] = 0.0;
  }

  while (iterations < maxIterations) {
    iterations++;
    //  The convergence threshold is based on the size of the elements in
    //  the strict upper triangle of the matrix.
    thresh = 0.0;
    for (int j = 0; j < n; j++) {
      for (int i = j + 1; i < n; i++) {
        thresh = thresh + pow(a(i, j), 2);
      }
    }

    thresh = sqrt(thresh) / static_cast<double>(4 * n);

    if (thresh == 0.0) {
      break;
    }

    for (int p = 0; p < n; p++) {
      for (int q = p + 1; q < n; q++) {
        gapq = 10.0 * fabs(a[p + q * n]);
        termp = gapq + fabs(eigenvals[p]);
        termq = gapq + fabs(eigenvals[q]);
        //  Annihilate tiny offdiagonal elements.
        if (4 < iterations && termp == fabs(eigenvals[p]) &&
            termq == fabs(eigenvals[p])) {
          a[p + q * n] = 0.0; // SEGMENTATION FAULT HERE
        }
        //  Otherwise, apply a rotation.
        else if (thresh <= fabs(a[p + q * n])) {
          h = eigenvals[q] - eigenvals[p];
          term = fabs(h) + gapq;

          if (term == fabs(h)) {
            t = a[p + q * n] / h;
          } else {
            theta = 0.5 * h / a[p + q * n];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) {
              t = -t;
            }
          }
          c = 1.0 / sqrt(1.0 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[p + q * n];
          //  Accumulate corrections to diagonal elements.
          zw[p] = zw[p] - h;
          zw[q] = zw[q] + h;
          eigenvals[p] -= h;
          eigenvals[q] += h;

          a[p + q * n] = 0.0;
          //  Rotate, using information from the upper triangle of A only.
          for (int j = 0; j < p; j++) {
            g = a[j + p * n];
            h = a[j + q * n];
            a[j + p * n] = g - s * (h + g * tau);
            a[j + q * n] = h + s * (g - h * tau);
          }

          for (int j = p + 1; j < q; j++) {
            g = a[p + j * n];
            h = a[j + q * n];
            a[p + j * n] = g - s * (h + g * tau);
            a[j + q * n] = h + s * (g - h * tau);
          }

          for (int j = q + 1; j < n; j++) {
            g = a[p + j * n];
            h = a[q + j * n];
            a[p + j * n] = g - s * (h + g * tau);
            a[q + j * n] = h + s * (g - h * tau);
          }
          //  Accumulate information in the eigenvector matrix.
          for (int j = 0; j < n; j++) {
            g = eigenvecs[j + p * n];
            h = eigenvecs[j + q * n];
            eigenvecs[j + p * n] = g - s * (h + g * tau);
            eigenvecs[j + q * n] = h + s * (g - h * tau);
          }
          rotations = rotations + 1;
        }
      }
    }

    for (int i = 0; i < n; i++) {
      bw[i] = bw[i] + zw[i];
      eigenvals[i] = bw[i];
      zw[i] = 0.0;
    }
  }
  //  Restore upper triangle of input matrix.
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++) {
      a[i + j * n] = a[j + i * n];
    }
  }
  //  Ascending sort the eigenvalues and eigenvectors.
  for (int k = 0; k < n - 1; k++) {
    int m = k;
    for (int l = k + 1; l < n; l++) {
      if (eigenvals[l] < eigenvals[m]) {
        m = l;
      }
    }

    if (m != k) {
      t = eigenvals[m];
      eigenvals[m] = eigenvals[k];
      eigenvals[k] = t;
      for (int i = 0; i < n; i++) {
        w = eigenvecs[i + m * n];
        eigenvecs[i + m * n] = eigenvecs[i + k * n];
        eigenvecs[i + k * n] = w;
      }
    }
  }

  delete[] bw;
  delete[] zw;

  return std::tuple<Matrix, Matrix, int, int>(eigenvecs, eigenvals, iterations,
                                              rotations);
}

// functions on VECTORS
double Matrix::dotProduct(Matrix a, Matrix b) {
  double sum = 0;
  for (int i = 0; i < a.rows_; ++i) {
    sum += (a(i, 0) * b(i, 0));
  }
  return sum;
}

// functions on AUGMENTED matrices
Matrix Matrix::augment(Matrix A, Matrix B) {
  Matrix AB(A.rows_, A.cols_ + B.cols_);
  for (int i = 0; i < AB.rows_; ++i) {
    for (int j = 0; j < AB.cols_; ++j) {
      if (j < A.cols_)
        AB(i, j) = A(i, j);
      else
        AB(i, j) = B(i, j - B.cols_);
    }
  }
  return AB;
}

Matrix Matrix::gaussianEliminate() {
  Matrix Ab(*this);
  int rows = Ab.rows_;
  int cols = Ab.cols_;
  int Acols = cols - 1;

  int i = 0; // row tracker
  int j = 0; // column tracker

  // iterate through the rows
  while (i < rows) {
    // find a pivot for the row
    bool pivot_found = false;
    while (j < Acols && !pivot_found) {
      if (Ab(i, j) != 0) { // pivot not equal to 0
        pivot_found = true;
      } else { // check for a possible swap
        int max_row = i;
        double max_val = 0;
        for (int k = i + 1; k < rows; ++k) {
          double cur_abs = Ab(k, j) >= 0 ? Ab(k, j) : -1 * Ab(k, j);
          if (cur_abs > max_val) {
            max_row = k;
            max_val = cur_abs;
          }
        }
        if (max_row != i) {
          Ab.swapRows(max_row, i);
          pivot_found = true;
        } else {
          j++;
        }
      }
    }

    // perform elimination as normal if pivot was found
    if (pivot_found) {
      for (int t = i + 1; t < rows; ++t) {
        for (int s = j + 1; s < cols; ++s) {
          Ab(t, s) = Ab(t, s) - Ab(i, s) * (Ab(t, j) / Ab(i, j));
          if (Ab(t, s) < EPS && Ab(t, s) > -1 * EPS)
            Ab(t, s) = 0;
        }
        Ab(t, j) = 0;
      }
    }

    i++;
    j++;
  }

  return Ab;
}

Matrix Matrix::rowReduceFromGaussian() {
  Matrix R(*this);
  int rows = R.rows_;
  int cols = R.cols_;

  int i = rows - 1; // row tracker
  int j = cols - 2; // column tracker

  // iterate through every row
  while (i >= 0) {
    // find the pivot column
    int k = j - 1;
    while (k >= 0) {
      if (R(i, k) != 0)
        j = k;
      k--;
    }

    // zero out elements above pivots if pivot not 0
    if (R(i, j) != 0) {

      for (int t = i - 1; t >= 0; --t) {
        for (int s = 0; s < cols; ++s) {
          if (s != j) {
            R(t, s) = R(t, s) - R(i, s) * (R(t, j) / R(i, j));
            if (R(t, s) < EPS && R(t, s) > -1 * EPS)
              R(t, s) = 0;
          }
        }
        R(t, j) = 0;
      }

      // divide row by pivot
      for (int k = j + 1; k < cols; ++k) {
        R(i, k) = R(i, k) / R(i, j);
        if (R(i, k) < EPS && R(i, k) > -1 * EPS)
          R(i, k) = 0;
      }
      R(i, j) = 1;
    }

    i--;
    j--;
  }

  return R;
}

/*void Matrix::readSolutionsFromRREF(ostream &os) {
  Matrix R(*this);

  // print number of solutions
  bool hasSolutions = true;
  bool doneSearching = false;
  int i = 0;
  while (!doneSearching && i < rows_) {
    bool allZeros = true;
    for (int j = 0; j < cols_ - 1; ++j) {
      if (R(i, j) != 0)
        allZeros = false;
    }
    if (allZeros && R(i, cols_ - 1) != 0) {
      hasSolutions = false;
      os << "NO SOLUTIONS" << endl << endl;
      doneSearching = true;
    } else if (allZeros && R(i, cols_ - 1) == 0) {
      os << "INFINITE SOLUTIONS" << endl << endl;
      doneSearching = true;
    } else if (rows_ < cols_ - 1) {
      os << "INFINITE SOLUTIONS" << endl << endl;
      doneSearching = true;
    }
    i++;
  }
  if (!doneSearching)
    os << "UNIQUE SOLUTION" << endl << endl;

  // get solutions if they exist
  if (hasSolutions) {
    Matrix particular(cols_ - 1, 1);
    Matrix special(cols_ - 1, 1);

    for (int i = 0; i < rows_; ++i) {
      bool pivotFound = false;
      bool specialCreated = false;
      for (int j = 0; j < cols_ - 1; ++j) {
        if (R(i, j) != 0) {
          // if pivot variable, add b to particular
          if (!pivotFound) {
            pivotFound = true;
            particular(j, 0) = R(i, cols_ - 1);
          } else { // otherwise, add to special solution
            if (!specialCreated) {
              special = Matrix(cols_ - 1, 1);
              specialCreated = true;
            }
            special(j, 0) = -1 * R(i, j);
          }
        }
      }
      os << "Special solution:" << endl << special << endl;
    }
    os << "Particular solution:" << endl << particular << endl;
  }
}*/

Matrix Matrix::inverse() {
  Matrix I = Matrix::createIdentity(rows_);
  Matrix AI = Matrix::augment(*this, I);
  Matrix U = AI.gaussianEliminate();
  Matrix IAInverse = U.rowReduceFromGaussian();
  Matrix AInverse(rows_, cols_);
  for (int i = 0; i < AInverse.rows_; ++i) {
    for (int j = 0; j < AInverse.cols_; ++j) {
      AInverse(i, j) = IAInverse(i, j + cols_);
    }
  }
  return AInverse;
}

Matrix Matrix::getMainDiagonal() {
  int minDim = std::min(rows_, cols_);
  Matrix res(1, minDim);
  for (int i = 0; i < minDim; i++) {
    res[i] = this->p[i][i];
  }
  return res;
}

/* PRIVATE HELPER FUNCTIONS
 ********************************/

void Matrix::allocSpace() {
  p = new double *[rows_];
  for (int i = 0; i < rows_; ++i) {
    p[i] = new double[cols_];
  }
}

Matrix Matrix::expHelper(const Matrix &m, int num) {
  if (num == 0) {
    return createIdentity(m.rows_);
  } else if (num == 1) {
    return m;
  } else if (num % 2 == 0) { // num is even
    return expHelper(m * m, num / 2);
  } else { // num is odd
    return m * expHelper(m * m, (num - 1) / 2);
  }
}

/* NON-MEMBER FUNCTIONS
 ********************************/

Matrix operator+(const Matrix &m1, const Matrix &m2) {
  Matrix temp(m1);
  return (temp += m2);
}

Matrix operator-(const Matrix &m1, const Matrix &m2) {
  Matrix temp(m1);
  return (temp -= m2);
}

Matrix operator*(const Matrix &m1, const Matrix &m2) {
  Matrix temp(m1);
  return (temp *= m2);
}

Matrix operator*(const Matrix &m, double num) {
  Matrix temp(m);
  return (temp *= num);
}

Matrix operator*(double num, const Matrix &m) { return (m * num); }

Matrix operator/(const Matrix &m, double num) {
  Matrix temp(m);
  return (temp /= num);
}

ostream &operator<<(ostream &os, const Matrix &m) {
  for (int i = 0; i < m.rows_; ++i) {
    os << m.p[i][0];
    for (int j = 1; j < m.cols_; ++j) {
      os << " " << m.p[i][j];
    }
    os << endl;
  }
  return os;
}

istream &operator>>(istream &is, Matrix &m) {
  for (int i = 0; i < m.rows_; ++i) {
    for (int j = 0; j < m.cols_; ++j) {
      is >> m.p[i][j];
    }
  }
  return is;
}