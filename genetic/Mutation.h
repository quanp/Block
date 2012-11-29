#ifndef GENETIC_MUTATION_H
#define GENETIC_MUTATION_H

#include <vector>
#include <algorithm>
#include "RandomGenerator.h"
#include "Gene.h"

namespace genetic
{

std::vector<int> PMUTATE(const std::vector<int>& base)
{
  const size_t& n = base.size();
  std::vector<int> mute(base);
  int ntry = irand(3) + 1;
  for(int k = 0; k < ntry; ++k) {
    int i = irand(n);
    int j = irand(n);
    if(i != j) std::swap(mute[i], mute[j]);
  }
  return mute;
}

std::vector<int> GMUTATE(const std::vector<int>& base)
{
  const size_t& n = base.size();
  std::vector<int> mute(n, 0);
  std::vector<int> ndiv(4, 0);
  for(int i = 0; i < 4; ++i) ndiv[i] = irand(n + 1);
  sort(ndiv.begin(), ndiv.end());
  int j = 0;
  for(int i = 0;       i < ndiv[0]; ++i) mute[j++] = base[i];
  for(int i = ndiv[2]; i < ndiv[3]; ++i) mute[j++] = base[i];
  for(int i = ndiv[1]; i < ndiv[2]; ++i) mute[j++] = base[i];
  for(int i = ndiv[0]; i < ndiv[1]; ++i) mute[j++] = base[i];
  for(int i = ndiv[3]; i < n;       ++i) mute[j++] = base[i];
  return mute;
}

//
// interfaces to class Gene
//
Gene PointMutate(const Gene& g)
{
  return Gene(PMUTATE(g.sequence()));
}

Gene GlobalMutate(const Gene& g)
{
  return Gene(GMUTATE(g.sequence()));
}

}; // namespace genetic

#endif // GENETIC_MUTATION_H
