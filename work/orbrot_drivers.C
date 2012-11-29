#include <cmath>
#include "orbrot_drivers.h"

void compute_onepdm  (const int& norbs,
                      const int& nelec,
                      const btas::DTensor<4>& twopdm,
                            btas::DTensor<2>& onepdm)
{
  double scale = 2.0 / (nelec - 1);
  if(!onepdm.data()) {
    onepdm.resize(norbs, norbs);
  }
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j < norbs; ++j) {
      onepdm(i, j) = 0.0;
      for(int k = 0; k < norbs; ++k) {
        onepdm(i, j) += twopdm(i, k, k, j);
      }
      onepdm(i, j) *= scale;
    }
  }
}

void compute_fock    (const int& norbs,
                      const btas::DTensor<2>& oneint,
                      const btas::DTensor<4>& twoint,
                      const btas::DTensor<2>& onepdm,
                      const btas::DTensor<4>& twopdm,
                            btas::DTensor<2>& fock)
{
  if(!fock.data()) {
    fock.resize(norbs, norbs);
  }
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j < norbs; ++j) {
      double sum = 0.0;
      for(int p = 0; p < norbs; ++p) {
        sum += oneint(i, p) * onepdm(j, p);
        for(int q = 0; q < norbs; ++q) {
          for(int r = 0; r < norbs; ++r) {
            sum += 2.0 * twoint(i, p, r, q) * twopdm(j, p, q, r);
          }
        }
      }
      fock(i, j) = sum;
    }
  }
}

void compute_gradient(const int& norbs,
                      const btas::DTensor<2>& fock,
                            btas::DTensor<2>& grad)
{
  if(!grad.data()) {
    grad.resize(norbs, norbs);
  }
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j < norbs; ++j) {
      grad(i, j) = 2.0 * (fock(i, j) - fock(j, i));
    }
  }
}

void compute_hessian (const int& norbs,
                      const btas::DTensor<2>& fock,
                      const btas::DTensor<2>& oneint,
                      const btas::DTensor<4>& twoint,
                      const btas::DTensor<2>& onepdm,
                      const btas::DTensor<4>& twopdm,
                            btas::DTensor<4>& hess)
{
  hess.resize(norbs, norbs, norbs, norbs);
  hess = 0.0;

  btas::DTensor<4> hw(norbs, norbs, norbs, norbs);
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j < norbs; ++j) {
      for(int k = 0; k < norbs; ++k) {
        hw(i, k, k, j) = fock(i, j) + fock(j, i);
      }
    }
  }
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j < norbs; ++j) {
      for(int k = 0; k < norbs; ++k) {
        for(int l = 0; l < norbs; ++l) {
          double wijkl = 0.0;
          for(int p = 0; p < norbs; ++p) {
            for(int q = 0; q < norbs; ++q) {
              wijkl += twopdm(k, p, q, i) * twoint(j, p, l, q);
              wijkl += twopdm(k, p, i, q) * twoint(j, l, p, q);
              wijkl += twopdm(k, i, p, q) * twoint(j, l, p, q);
            }
          }
          hw(i, j, k, l) += 2.0 * onepdm(i, k) * oneint(j, l) + 4.0 * wijkl;
        }
      }
    }
  }
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j < norbs; ++j) {
      for(int k = 0; k < norbs; ++k) {
        for(int l = 0; l < norbs; ++l) {
          hess(i, j, k, l) = hw(i, j, k, l) - hw(j, i, k, l) - hw(i, j, l, k) + hw(j, i, l, k);
        }
      }
    }
  }
}

double compute_energy(const int& norbs,
                      const btas::DTensor<2>& oneint,
                      const btas::DTensor<4>& twoint,
                      const btas::DTensor<2>& onepdm,
                      const btas::DTensor<4>& twopdm)
{
  double energy = 0.0;
  for(int p = 0; p < norbs; ++p) {
    for(int q = 0; q < norbs; ++q) {
      energy += oneint(p, q) * onepdm(p, q);
      for(int r = 0; r < norbs; ++r) {
        for(int s = 0; s < norbs; ++s) {
          energy += twoint(p, q, s, r) * twopdm(p, q, r, s);
        }
      }
    }
  }
  return energy;
}

void compute_rotation(const btas::DTensor<2>& x, btas::DTensor<2>& u)
{
  btas::DTensor<2> x2;
  btas::Dgemm(btas::NoTrans, btas::NoTrans, 1.0, x, x, 1.0, x2);
  btas::Dscal(-1.0, x2);
  btas::DTensor<1> t;
  btas::DTensor<2> p;
  btas::Dsyev(x2, t, p);
  for(int i = 0; i < t.size(); ++i) {
    if(t(i) > 1.0e-16) {
      t(i)   = sqrt(t(i));
      t(i+1) = t(i);
      ++i;
    }
    else {
      t(i)   = 0.0;
    }
  }
  btas::DTensor<1> cos_t(t.shape());
  btas::DTensor<1> sin_t(t.shape());
  for(int i = 0; i < t.size(); ++i) {
    if(t(i) > 1.0e-8) {
      cos_t(i) = cos(t(i));
      sin_t(i) = sin(t(i)) / t(i);
    }
    else {
      cos_t(i) = 1.0;
      sin_t(i) = 0.0;
    }
  }

  u.free();

  btas::DTensor<2> cp;
  btas::Dcopy(p, cp);
  btas::Dright_update(cos_t, cp);
  btas::Dgemm(btas::Trans, btas::NoTrans, 1.0, p, cp, 1.0, u);

  btas::DTensor<2> sp;
  btas::Dcopy(p, sp);
  btas::Dright_update(sin_t, sp);
  btas::DTensor<2> psp;
  btas::Dgemm(btas::Trans,   btas::NoTrans, 1.0, p, sp, 1.0, psp);
  btas::Dgemm(btas::NoTrans, btas::NoTrans, 1.0, psp, x, 1.0, u);
}
