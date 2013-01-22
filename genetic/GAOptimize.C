#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "Evaluate.h"
#include "Generation.h"
#include "GAInput.h"
#include "SymmetricMatrix.h"
#include "DistanceMatrix.h"
using namespace std;

#include <boost/function.hpp>
#include <boost/bind.hpp>
using boost::function;
using boost::bind;

#ifndef SERIAL
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

genetic::GAInput gainput;

genetic::Cell comm_optimize(void)
{
  using namespace genetic;
  Generation ancestor;
//fout << "--------------------------- STARTING GA OPTIMIZATION ----------------------------" << endl;
  for(int g = 0; g < gainput.max_generation; ++g) {
    Generation nextgen;
    nextgen.generate(ancestor);
    ancestor = nextgen;
//  if(g % 1000 == 0) fout << "Generation [ " << setw(8) << g << " ]: " << ancestor.min() << endl;
  }
//fout << "-------------------- RESULTS: POPULATION OF FINAL GENERATION --------------------" << endl;
//fout << ancestor << endl;
  return ancestor.min();
}

genetic::Cell GAOptimize(std::ifstream& confFile, std::ifstream& dumpFile)
{
  using namespace genetic;
#ifndef SERIAL
  mpi::communicator world;
#endif
  DistanceMatrix Dist;
  SymmetricMatrix Kint;
  double Ksum = 0.0;
#ifndef SERIAL
  if(world.rank() == 0)
  {
#endif
    if(confFile.is_open()) gainput.configure(confFile);

    Kint = read_integral(dumpFile);
    Gene::Length() = Kint.size();

    if(gainput.graph.size() > 0) {
      if(gainput.graph.size() != Gene::Length()) {
        cout << "number of sites in config file mismatched to that in integral file" << endl;
        abort();
      }
      Dist.reset(gainput.scale, gainput.graph);
    }
    else {
      std::vector<int> mps_graph(Gene::Length(),-1);
      for(int i = 0; i < Gene::Length(); ++i) mps_graph[i] = i - 1;
      Dist.reset(gainput.scale, mps_graph);
    }

    if(gainput.max_cells == 0) gainput.max_cells = 2 * Gene::Length();
    for(int i = 0; i < Kint.size(); ++i)
      for(int j = 0; j < i; ++j) Ksum += Kint(i, j);

#ifndef SERIAL
  }
  mpi::broadcast(world, gainput, 0);
  mpi::broadcast(world, Gene::Length(), 0);
  mpi::broadcast(world, Dist, 0);
  mpi::broadcast(world, Kint, 0);
  mpi::broadcast(world, Ksum, 0);
#endif
  Cell::cost_functor = bind(Evaluate, 1.0/Ksum, _1, Dist, Kint);

  Cell best;
#ifndef SERIAL
  int nproc = world.size();
  int nrank = world.rank();
  int ntask = 1 + gainput.max_community / nproc;
  srand(gainput.random_seed + nrank);

  Cell comm_best = comm_optimize();
  cout << "Order #" << setw(3) << nrank << ": " << comm_best << endl;
  for(int i = 1; i < ntask; ++i) {
    Cell comm_cell = comm_optimize();
    cout << "Order #" << setw(3) << i * nproc + nrank << ": " << comm_cell << endl;
    if(comm_cell < comm_best) comm_best = comm_cell;
  }
  if(nrank == 0)
    mpi::reduce(world, comm_best, best, mpi::minimum<Cell>(), 0);
  else
    mpi::reduce(world, comm_best,       mpi::minimum<Cell>(), 0);
#else
  int ntask = gainput.max_community;
  srand(gainput.random_seed);

  best = comm_optimize();
  cout << "Order #" << setw(3) << 0 << ": " << best << endl;
  for(int i = 1; i < ntask; ++i) {
    Cell comm_cell = comm_optimize();
    cout << "Order #" << setw(3) << i << ": " << comm_cell << endl;
    if(comm_cell < best) best = comm_cell;
  }
#endif

  return best;
}

std::vector<int> gaorder(std::ifstream& confFile, std::ifstream& dumpFile)
{
#ifndef SERIAL
  mpi::communicator world;
#endif
  genetic::Cell final = GAOptimize(confFile, dumpFile);
#ifndef SERIAL
  if(world.rank() == 0)
#endif
  {
    cout << "##################### MINIMUM GENE REP. #####################" << endl;
    cout << "Gene with MinValue = " << final << endl;
    cout << "Effective Distance = " << sqrt(final.fit()) << endl;
  }

  return final.gene().sequence();
}
