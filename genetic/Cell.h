#ifndef GENETIC_CELL_H
#define GENETIC_CELL_H

#include <iostream>
#include <boost/function.hpp>
#include <boost/serialization/serialization.hpp>
#include "Gene.h"

namespace genetic
{

//
// class Cell defines individuals in genetic algorithm
//
class Cell
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & m_gen;
    ar & m_fit;
  }
private:
  Gene
    m_gen;
  double
    m_fit;
public:
  static boost::function<double(const std::vector<int>&)> cost_functor;

  Cell()
  {
  }
  Cell(const Gene& gen)
  {
    create(gen);
  }
  void create(const Gene& gen)
  {
    m_gen = gen;
    m_fit = cost_functor(gen.sequence());
  }

  inline bool operator== (const Cell& other) const { return m_fit == other.m_fit; }
  inline bool operator!= (const Cell& other) const { return m_fit != other.m_fit; }
  inline bool operator<  (const Cell& other) const { return m_fit <  other.m_fit; }
  inline bool operator>  (const Cell& other) const { return m_fit >  other.m_fit; }

  const Gene& gene() const
  {
    return m_gen;
  }
  const double& fit() const
  {
    return m_fit;
  }
  friend std::ostream& operator<< (std::ostream& ost, const Cell& cell)
  {
    ost.setf(std::ios::fixed, std::ios::floatfield);
    ost.precision(2);
    ost << cell.m_gen << " : fit = " << cell.m_fit;
    return ost;
  }
};

boost::function<double(const std::vector<int>&)> Cell::cost_functor;

}; // namespace genetic

#endif // GENETIC_CELL_H
