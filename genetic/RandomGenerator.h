#ifndef GENETIC_RANDOM_GENERATOR_H
#define GENETIC_RANDOM_GENERATOR_H

#include <vector>
#include <map>
#include <cstdlib>

namespace genetic
{

// random number generator
inline int irand(const int& n)
{
  int r = (int) (((double)(rand())) * n / RAND_MAX);
  return (r < n) ? r : n - 1;
}

inline int irand(const int& nmin, const int& nmax)
{
  int n = nmax - nmin;
  return nmin + irand(n);
}

inline double drand(const double& alpha) { return (alpha * (((double)(rand())) / RAND_MAX)); }

// random bit
inline int brand(const double& alpha)
{
  if(drand(1.0) < alpha) return 1;
  else                   return 0;
}

// Random Bit Sequence Generator
std::vector<int> RandomBitSequence(const double& alpha, const int& n)
{
  std::vector<int> seq(n, 0);
  for(int i = 0; i < n; ++i) seq[i] = brand(alpha);
  return seq;
}

std::vector<int> RandomBitSequence(const int& n)
{
  return RandomBitSequence(0.5, n);
}

// Random Sequence Generator
std::vector<int> RandomSequence(const int& n)
{
  std::vector<int> seq(n, 0);
  std::multimap<double, int> rand_map;
  for(int i = 0; i < n; ++i) rand_map.insert(std::make_pair(drand(1.0), i));
  std::multimap<double, int>::iterator it = rand_map.begin();
  for(int i = 0; i < n; ++i, ++it) seq[i] = it->second;
  return seq;
}

}; // namespace genetic

#endif // GENETIC_RANDOM_GENERATOR_H
