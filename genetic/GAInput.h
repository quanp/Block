#ifndef GENETIC_GAINPUT_H
#define GENETIC_GAINPUT_H

#include <fstream>
#include <cstring>
#include <vector>
#include <ctime>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/algorithm/string.hpp>
#include "SymmetricMatrix.h"

namespace genetic
{

enum SELECTION_TYPE { GAUSS, BOLTZMAN, UNIFORM };

//
// parsing molpro FCIDUMP file to return integral object
//

std::vector<std::string> gettoken(std::ifstream& fin)
{
  // read msg from fin
  std::string msg;
  std::getline(fin, msg);
  boost::trim(msg);
  // split msg into token
  std::vector<std::string> tok;
  boost::split(tok, msg, boost::is_any_of("=, \t"), boost::token_compress_on);
  return tok;
}

SymmetricMatrix read_integral(std::ifstream& fdump)
{
  SymmetricMatrix moint;
  std::vector<std::string> tok;

  int norbs = 0;
  // read first line
  tok = gettoken(fdump);
  for(int i = 0; i < tok.size(); ++i) {
    if(tok[i] == "NORB" ) norbs = atoi(tok[++i].c_str());
  }
  moint.resize(norbs);
  // find &END control
  while(1) {
    tok = gettoken(fdump);
    if(tok[0] == "&END" || tok[0] == "/") break;
  }
  // read integrals
  while((tok = gettoken(fdump)).size() > 1) {
    double value = atof(tok[0].c_str());
    int i = atoi(tok[1].c_str()) - 1;
    int j = atoi(tok[2].c_str()) - 1;
    int k = atoi(tok[3].c_str()) - 1;
    int l = atoi(tok[4].c_str()) - 1;
    // exchange integrals
    if(i > 0 && j > 0 && i == k && j == l) {
      moint(i, j) = value;
    }
  }
  return moint;
}

//
// input class
//
struct GAInput
{
private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & graph;
    ar & max_community;
    ar & max_generation;
    ar & max_cells;
    ar & thre_cloning;
    ar & thre_mutation;
    ar & max_elite;
    ar & scale;
    ar & weight;
    ar & random_seed;
    ar & select;
  }
public:
  std::vector<int> graph;
  int              max_community;
  int              max_generation;
  int              max_cells;
  double           thre_cloning;
  double           thre_mutation;
  int              max_elite;
  double           scale;
  double           weight;
  unsigned int     random_seed;
  SELECTION_TYPE   select;

  GAInput()
  {
    max_community  = 10;
    max_generation = 10000;
    max_cells      = 0;
    thre_cloning   = 0.90;
    thre_mutation  = 0.10;
    max_elite      = 1;
    scale          = 1.0;
    weight         = 1.0;
    select         = GAUSS;
    random_seed    = time(NULL);
  }

  GAInput(std::ifstream& config)
  {
    max_community  = 10;
    max_generation = 10000;
    max_cells      = 0;
    thre_cloning   = 0.90;
    thre_mutation  = 0.10;
    max_elite      = 1;
    scale          = 1.0;
    weight         = 1.0;
    select         = GAUSS;
    random_seed    = time(NULL);

    configure(config);
  }

  void configure(std::ifstream& config)
  {
    graph.clear();
    std::string entry;
    while(config >> entry) {
      if(entry == "graph") {
        int n_vert;
        config >> n_vert;
        graph = std::vector<int>(n_vert,-1);
        for(int i = 0; i < n_vert; ++i) {
          int iconn;
          config >> iconn;
          graph[i] = iconn - 1;
        }
      }
      if(entry == "maxcomm")  config >> max_community;
      if(entry == "maxgen")   config >> max_generation;
      if(entry == "maxcell")  config >> max_cells;
      if(entry == "cloning")  config >> thre_cloning;
      if(entry == "mutation") config >> thre_mutation;
      if(entry == "elite")    config >> max_elite;
      if(entry == "scale")    config >> scale;
      if(entry == "weight")   config >> weight;
      if(entry == "seed")     config >> random_seed;
      if(entry == "select" || entry == "method") {
        config >> entry;
        if(entry == "gauss")    select = GAUSS;
        if(entry == "boltzman") select = BOLTZMAN;
        if(entry == "uniform")  select = UNIFORM;
      }
    }
  }
};

};

#endif // GENETIC_GAINPUT_H
