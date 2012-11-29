#ifndef GENETIC_GENERATION_H
#define GENETIC_GENERATION_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include "Cell.h"
#include "GAInput.h"

extern genetic::GAInput gainput;

namespace genetic
{

//
// generation
//
class Generation
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & m_cells;
    ar & m_prob;
  }

  void ComputeProbability()
  {
    // sort by fitness
    sort(m_cells.begin(), m_cells.end());
         if(gainput.select == GAUSS)    GaussProb();
    else if(gainput.select == BOLTZMAN) BoltzmanProb();
    else                                UniformProb();
  }

  void GaussProb()
  {
    const int& n = m_cells.size();
    m_prob = std::vector<double>(n, 0.0);
    for(int i = 0; i < n; ++i) {
      m_prob[i] = sqrt(m_cells[i].fit());
    }
    double minfit = m_prob[0];
    for(int i = 0; i < n; ++i) {
      double value = m_prob[i] - minfit;
      m_prob[i] = exp(-gainput.weight * value * value);
    }
  }
  void BoltzmanProb(void)
  {
    const int& n = m_cells.size();
    m_prob = std::vector<double>(n, 0.0);
    for(int i = 0; i < n; ++i) {
      m_prob[i] = sqrt(m_cells[i].fit());
    }
    double minfit = m_prob[0];
    for(int i = 0; i < n; ++i) {
      double value = m_prob[i] - minfit;
      m_prob[i] = exp(-gainput.weight * value);
    }
  }
  void UniformProb(void)
  {
    const int& n = m_cells.size();
    m_prob = std::vector<double>(n, 0.0);
    for(int i = 0; i < n; ++i) {
      m_prob[i] = 1.0;
    }
  }
public:
  Generation()
  {
    initialize();
  }
  // God makes life
  void initialize()
  {
    const int& n = gainput.max_cells;
    m_cells = std::vector<Cell>(n, Cell());
    for(int i = 0; i < n; ++i) {
      m_cells[i].create(Gene::Random());
    }
    ComputeProbability();
  }

  const Cell& select() const
  {
    const int& n = m_cells.size();
    int i = 0;
    while(1) {
      i = irand(n);
      if(m_prob[i] > drand(1.0)) break;
    }
    return m_cells[i];
  }
  Cell& select()
  {
    const int& n = m_cells.size();
    int i = 0;
    while(1) {
      i = irand(n);
      if(m_prob[i] > drand(1.0)) break;
    }
    return m_cells[i];
  }

  inline const Cell& select(const int& i) const { return m_cells[i]; }
  inline       Cell& select(const int& i)       { return m_cells[i]; }

  inline const Cell& min() const { return *min_element(m_cells.begin(), m_cells.end()); }
  inline       Cell& min()       { return *min_element(m_cells.begin(), m_cells.end()); }

  inline const Cell& max() const { return *max_element(m_cells.begin(), m_cells.end()); }
  inline       Cell& max()       { return *max_element(m_cells.begin(), m_cells.end()); }

  void generate(const Generation& ancestor)
  {
    const int& n = ancestor.m_cells.size();
    m_cells = std::vector<Cell>(n, Cell());

    int i = 0;
    // chose elites
    for(; i < gainput.max_elite; ++i) m_cells[i] = ancestor.m_cells[i];
    // selection to next generation
    while(i < n) {
      Gene child;
      if(drand(1.0) > gainput.thre_cloning) {
        child = ancestor.select().gene();
      }
      else {
        Cell father(ancestor.select());
        Cell mother(ancestor.select());
        child = father.gene() * mother.gene();
      }
      // mutation
      if(drand(1.0) < gainput.thre_mutation) child.p_mutate();
      if(drand(1.0) < gainput.thre_mutation) child.g_mutate();

      m_cells[i++] = Cell(child);
    }
    ComputeProbability();
  }

  friend std::ostream& operator<< (std::ostream& ost, const Generation& g)
  {
    using std::setw;
    using std::endl;
    ost.setf(std::ios::fixed, std::ios::floatfield);
    ost.precision(4);
    for(int i = 0; i < g.m_cells.size(); ++i) {
      ost << setw(4) << i << ": " << g.m_cells[i] << " (" << g.m_prob[i] << ")" << endl;
    }
    return ost;
  }

private:
  std::vector<Cell>
    m_cells;
  std::vector<double>
    m_prob;
};

}; // namespace genetic

#endif // GENETIC_GENERATION_H
