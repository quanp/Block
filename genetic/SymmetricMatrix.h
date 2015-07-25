#ifndef GENETIC_SYMMETRIC_MATRIX_H
#define GENETIC_SYMMETRIC_MATRIX_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace genetic
{

//
// symmetric matrix
//
class SymmetricMatrix
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & m_size;
    ar & m_data;
  }
protected:
  // indexing function of symmetric matrix
  inline int Indexing(const int& i, const int& j) const
  {
    return i < j ? ((j * (j + 1)) / 2 + i):((i * (i + 1)) / 2+ j);
  }
public:
  SymmetricMatrix()
  {
    m_size = 0;
  }
  SymmetricMatrix(const size_t& size)
  {
    resize(size);
  }
  void resize(const size_t& size)
  {
    m_size = size;
    m_data = std::vector<double>(size * (size + 1) / 2, 0.0);
  }
  inline size_t size() const
  {
    return m_size;
  }
  inline const double& operator() (const int& i, const int& j) const
  {
    return m_data[Indexing(i, j)];
  }
  inline double& operator() (const int& i, const int& j)
  {
    return m_data[Indexing(i, j)];
  }
  // print out matrix
  friend std::ostream& operator<< (std::ostream& ost, const SymmetricMatrix& data)
  {
    using std::setw;
    using std::endl;
    const int ncols = 10;
    ost.setf(std::ios::fixed, std::ios::floatfield);
    ost.precision(8);
    ost << "\tsize = " << setw(3) << data.size() << endl;

    int ndiv = data.size() / ncols;
    for(int t = 0; t < ndiv; ++t) {
      int n = t * ncols;
      for(int j = n; j < data.size(); ++j) {
        ost << "\t";
        for(int i = n; i <= std::min(j, n + ncols); ++i) ost << setw(12) << data(i, j);
        ost << endl;
      }
    }
    for(int j = ndiv * ncols; j < data.size(); ++j) {
      ost << "\t";
      for(int i = ndiv * ncols; i <= j; ++i) ost << setw(12) << data(i, j);
      ost << endl;
    }
    return ost;
  }
protected:
  size_t
    m_size;
  std::vector<double>
    m_data;
};

}; // namespace genetic

#endif // GENETIC_SYMMETRIC_MATRIX_H
