#ifndef ORBROT_AUGMENTED_HESSIAN_H
#define ORBROT_AUGMENTED_HESSIAN_H

//
// input : hessian (matrix), gradient (vector)
// output: displacement (vector)
//
void aughessian(const Matrix& hess, const ColumnVector& grad, ColumnVector& disp);

#endif // ORBROT_AUGMENTED_HESSIAN_H
