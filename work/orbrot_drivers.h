#ifndef ORBROT_DRIVERS_H
#define ORBROT_DRIVERS_H

#include "orbrot_defs.h"

void compute_onepdm  (const int& norbs,
                      const int& nelec,
                      const Matrix& twopdm,
                            ColumnVector& onepdm);

void compute_fock    (const int& norbs,
                      const ColumnVector& oneint,
                      const Matrix& twoint,
                      const ColumnVector& onepdm,
                      const Matrix& twopdm,
                            ColumnVector& fock);

void compute_gradient(const int& norbs,
                      const ColumnVector& fock,
                            ColumnVector& grad);

void compute_hessian (const int& norbs,
                      const ColumnVector& fock,
                      const ColumnVector& oneint,
                      const Matrix& twoint,
                      const ColumnVector& onepdm,
                      const Matrix& twopdm,
                            Matrix& hess);

double compute_energy(const int& norbs,
                      const ColumnVector& oneint,
                      const Matrix& twoint,
                      const ColumnVector& onepdm,
                      const Matrix& twopdm);

void compute_rotation(const ColumnVector& x,
                            ColumnVector& u);

#endif // ORBROT_DRIVERS_H
