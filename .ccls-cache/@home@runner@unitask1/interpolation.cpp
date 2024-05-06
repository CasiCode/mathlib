#include <cmath>

#include "interpolation.hpp"

#define PI 3.14159265358979323846

double *chebyshevCoefficients(double a, double b, int n, double f(double x)) {

  //  Purpose:
  //    CHEBYSHEV_COEFFICIENTS determines Chebyshev interpolation coefficients.
  //  Reference:
  //    Roger Broucke,
  //    Algorithm 446:
  //    Ten Subroutines for the Manipulation of Chebyshev Series,
  //    Communications of the ACM,
  //    Volume 16, Number 4, April 1973, pages 254-256.
  //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
  //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
  //    Second Edition,
  //    Cambridge University Press, 1992,
  //    ISBN: 0-521-43064-X,
  //    LC: QA297.N866.
  //  Parameters:
  //    Input, double A, B, the domain of definition.
  //    Input, int N, the order of the interpolant.
  //    Input, double F ( double X ), an external function.
  //    Output, double CHEBYSHEV_COEFFICIENTS[N], the Chebyshev coefficients.

  double angle;
  double *c;
  double *fx;
  int i;
  int j;
  double x;

  fx = new double[n];

  for (i = 0; i < n; i++) {
    angle = static_cast<double>(2 * i + 1) * PI / static_cast<double>(2 * n);
    x = cos(angle);
    x = 0.5 * (a + b) + x * 0.5 * (b - a);
    fx[i] = f(x);
  }

  c = new double[n];

  for (i = 0; i < n; i++) {
    c[i] = 0.0;
    for (j = 0; j < n; j++) {
      angle = static_cast<double>(i * (2 * j + 1)) * PI / static_cast<double>(2 * n);
      c[i] = c[i] + fx[j] * cos(angle);
    }
  }

  for (i = 0; i < n; i++) {
    c[i] = 2.0 * c[i] / static_cast<double>(n);
  }

  delete[] fx;

  return c;
}

double *chebyshevInterpolant(double a, double b, int n, double c[], int m,
                              double x[]) {
  //  Purpose:
  //    CHEBYSHEV_INTERPOLANT evaluates a Chebyshev interpolant.
  //  Reference:
  //    Roger Broucke,
  //    Algorithm 446:
  //    Ten Subroutines for the Manipulation of Chebyshev Series,
  //    Communications of the ACM,
  //    Volume 16, Number 4, April 1973, pages 254-256.
  //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
  //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
  //    Second Edition,
  //    Cambridge University Press, 1992,
  //    ISBN: 0-521-43064-X,
  //    LC: QA297.N866.
  //  Parameters:
  //    Input, double A, B, the domain of definition.
  //    Input, int N, the order of the polynomial.
  //    Input, double C[N], the Chebyshev coefficients.
  //    Input, int M, the number of points.
  //    Input, double X[M], the point at which the polynomial is
  //    to be evaluated.
  //    Output, double CHEBYSHEF_INTERPOLANT[M], the value of the Chebyshev
  //    polynomial at X.

  double *cf;
  double di;
  double dip1;
  double dip2;
  int i;
  int j;
  double y;

  cf = new double[m];

  for (j = 0; j < m; j++) {
    dip1 = 0.0;
    di = 0.0;
    y = (2.0 * x[j] - a - b) / (b - a);

    for (i = n - 1; 1 <= i; i--) {
      dip2 = dip1;
      dip1 = di;
      di = 2.0 * y * dip1 - dip2 + c[i];
    }
    cf[j] = y * di - dip1 + 0.5 * c[0];
  }

  return cf;
}

double *chebyshevZeros(int n) {
  //  Purpose:
  //    CHEBYSHEV_ZEROS returns zeroes of the Chebyshev polynomial T(N)(X).
  //  Discussion:
  //    We produce the Chebyshev zeros in ascending order.
  //  Parameters:
  //    Input, int N, the order of the polynomial.
  //    Output, double CHEBYSHEV_ZEROS[N], the zeroes of T(N)(X).
  double angle;
  int i;
  double *x;

  x = new double[n];

  for (i = 0; i < n; i++) {
    angle = static_cast<double>(2 * (n - i) - 1) * PI /
            static_cast<double>(2 * n);
    x[i] = cos(angle);
  }

  return x;
}

// function to interpolate the given data points using Lagrange's formula
// xi corresponds to the new data point whose value is to be obtained
double *lagrangeInterpolation(Nodes data, int n, double x[]) {
    double *res;

    res = new double[n];
    for (int k = 0; k < n; k++) { //might not be the correct placement of the for loop !!!
        double result = 0; // Initialize result
        for (int i = 0; i < data.size(); i++) {
            // Compute individual terms of above formula
            double term = data[i].y();
            for (int j = 0; j < n; j++) {
                if (j != i) { term = term * (x[k] - data[j].x()) / double(data[i].x() - data[j].x()); }
            }

            // Add current term to result
            result += term;
        }
        res[k] = result;
    }
    return res;
}

double *chebyshevInterpolation(double a, double b, int n, int m, double x[], double f(double x)) {
    double *coeffs = chebyshevCoefficients(a, b, n, f);
    double *res = chebyshevInterpolant(a, b, n, coeffs, m, x);
    delete[] coeffs;
    return res;
}