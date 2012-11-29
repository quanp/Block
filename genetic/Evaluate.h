#ifndef GENETIC_EVALUATE_H
#define GENETIC_EVALUATE_H

#include <vector>
#include "SymmetricMatrix.h"

namespace genetic
{

//
// permute integral by given ordering
//
void Permute(const SymmetricMatrix& tin, const std::vector<int>& sequence, SymmetricMatrix& tout)
{
  const size_t& n = sequence.size();
  tout.resize(n);
  for(int i = 0; i < n; ++i) {
    int ip = sequence[i];
    for(int j = 0; j <= i; ++j) {
      int jp = sequence[j];
      tout(i, j) = tin(ip, jp);
    }
  }
}

//
// definition of cost function
//
double Evaluate(const double& scale, const std::vector<int>& sequence, const SymmetricMatrix& dist, const SymmetricMatrix& t)
{
  const size_t& n = sequence.size();
  SymmetricMatrix tp;
  Permute(t, sequence, tp);
  double fcost = 0.0;
  for(int i = 1; i < n; ++i) {
    for(int j = 0; j < i; ++j) {
      int d = dist(i, j);
      fcost += d * d * tp(i, j);
    }
  }
  return scale * fcost;
}

}; // namespace genetic

#endif // GENETIC_EVALUATE_H
