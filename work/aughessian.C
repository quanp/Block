#include <iostream>
#include <iomanip>
#include "aughessian.h"

void aughessian(const SymmetricMatrix& hess, const ColumnVector& grad, ColumnVector& disp)
{
  int n = hess.Nrows();
  // construct aug-hessian matrix
  // [ 0  g ]
  // [ g  h ]
  SymmetricMatrix augh(n+1);
  augh(0, 0) = 0.0;
  for(int j = 0; j < n; ++j)
    augh(1+j, 0) = grad(j);
  for(int i = 0; i < n; ++i)
    for(int j = 0; j <= i; ++j)
      augh(1+i, 1+j) = hess(i, j);
  // diagonalize augh
  DiagonalMatrix eigs;
  Matrix         augw;
  SpinAdapted::diagonalise(augh, eigs, augw);

  disp.ReSize(n);
  disp = 0.0;
  bool stationary = true;
  for(int i = 0; i < n+1; ++i) {
    double a0 = augw(i, 0);
    if(eigs(i) < 0.0 && fabs(a0) > 1.0e-8) {
      stationary = false;
      for(int j = 0; j < n; ++j)
        disp(j) = augw(i, 1+j) / a0;
      break;
    }
  }
}
