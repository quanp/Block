#ifndef GENETIC_DISTANCE_MATRIX_H
#define GENETIC_DISTANCE_MATRIX_H

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "SymmetricMatrix.h"

namespace genetic
{

//
// distance matrix which represents weighted graph
//
class DistanceMatrix : public SymmetricMatrix
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<SymmetricMatrix>(*this);
  }
public:
  DistanceMatrix() : SymmetricMatrix()
  {
  }
  DistanceMatrix(const std::vector<int>& connect)
  {
    reset(connect);
  }
  DistanceMatrix(const double& scale, const std::vector<int>& connect)
  {
    reset(scale, connect);
  }
  inline size_t size() const
  {
    return m_size;
  }
  inline const double& operator() (const int& i, const int& j) const
  {
    return m_data[Indexing(i, j)];
  }
  void reset(const std::vector<int>& connect)
  {
    reset(1.0, connect);
  }
  //
  // create weighted graph
  //
  void reset(const double& scale, const std::vector<int>& connect)
  {
    SymmetricMatrix::resize(connect.size());
    // compute size of entanglement
    std::vector<int> ndims(m_size, 1);
    for(int i = m_size - 1; i > 0; --i) {
      ++ndims[connect[i]];
    }
    for(int i = 0; i < m_size; ++i) {
      if(ndims[i] > m_size - ndims[i]) ndims[i] = m_size = ndims[i];
    }
    // compute weighted distance matrix step 1: setting neighboring pair
    for(int i = 1; i < m_size; ++i) {
      double dist = 1.0;
      for(int j = 0; j < ndims[i]; ++j) dist *= scale;
      m_data[Indexing(i, connect[i])] = dist;
    }
    // compute weighted distance matrix step 2: compute all pair recursively
    for(int i = 1; i < m_size; ++i) {
      int k = connect[i];
      for(int j = 0; j < i; ++j) {
        if(j == k) continue;
        m_data[Indexing(i, j)] = m_data[Indexing(i, k)] + m_data[Indexing(j, k)];
      }
    }
  }
};

}; // namespace genetic

#endif // GENETIC_DISTANCE_MATRIX_H
