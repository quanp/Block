#ifndef GENETIC_GENE_H
#define GENETIC_GENE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "RandomGenerator.h"

namespace genetic
{

class Gene;
Gene CrossOver(const Gene&, const Gene&);
Gene PointMutate(const Gene&);
Gene GlobalMutate(const Gene&);

//
// gene expression class
//
class Gene
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & m_sequence;
  }
private:
  static int
    m_length;
  std::vector<int>
    m_sequence;
public:
  static int& Length()
  {
    return m_length;
  }
  const static Gene Random()
  {
    return Gene(RandomSequence(m_length));
  }
  // constructor
  Gene()
  {
  }
  Gene(const std::vector<int>& sequence)
  {
    if(m_length == 0) m_length = sequence.size();
    m_sequence = sequence;
  }
  const std::vector<int>& sequence() const
  {
    return m_sequence;
  }
  Gene  operator* (const Gene& other) const
  {
    return Gene(CrossOver(*this, other));
  }
  void p_mutate() // point mutation
  {
    *this = PointMutate(*this);
  }
  void g_mutate() // global mutation
  {
    *this = GlobalMutate(*this);
  }
  friend std::ostream& operator<< (std::ostream& ost, const Gene& g)
  {
    int n = m_length - 1;
    ost << "sequence = {";
    for(int i = 0; i < n; ++i)
    ost << std::setw(3) << g.m_sequence[i] + 1 << ",";
    ost << std::setw(3) << g.m_sequence[n] + 1 << "}";
    return ost;
  }
};

int Gene::m_length = 0;

}; // namespace genetic

#include "CrossOver.h"
#include "Mutation.h"

#endif // GENETIC_GENE_H
