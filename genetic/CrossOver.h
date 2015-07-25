#ifndef GENETIC_CROSS_OVER_H
#define GENETIC_CROSS_OVER_H

#include <vector>
#include "RandomGenerator.h"
#include "Gene.h"

namespace genetic
{

//
// ordered crossover ( Partically Mapped Crossover-Uniform, PMXU )
//
std::vector<int> PMXU(const std::vector<int>& seq1, const std::vector<int>& seq2)
{
  int n = seq1.size();
  std::vector<int> seq1_map(n, 0);
  std::vector<int> seq2_map(n, 0);
  for(int i = 0; i < n; ++i) {
    seq1_map[seq1[i]] = i;
    seq2_map[seq2[i]] = i;
  }
  std::vector<int> mask(RandomBitSequence(n));
  std::vector<int> seq3(n, -1);

  std::vector<int> tmp2(seq2);
  for(int i = 0; i < n; ++i) {
    if(mask[i] == 1) {
      seq3[i] = seq1[i];
      int k = seq1_map[seq2[i]];
      if(mask[k] == 0) {
        int j = seq2_map[seq1[i]];
        while(mask[j]) j = seq2_map[seq1[j]];
        std::swap(tmp2[i], tmp2[j]);
      }
    }
  }
  for(int i = 0; i < n; ++i) if(mask[i] == 0) seq3[i] = tmp2[i];
  return seq3;
}

//
// interface to class Gene
//
Gene CrossOver(const Gene& g1, const Gene& g2)
{
  return Gene(PMXU(g1.sequence(), g2.sequence()));
}

}; // namespace genetic

#endif // GENETIC_CROSS_OVER_H
