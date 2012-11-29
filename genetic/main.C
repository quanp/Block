#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
using namespace std;

#ifndef SERIAL
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

vector<int> gaorder(ifstream&, ifstream&);

int main(int argc, char* argv[])
{
#ifndef SERIAL
  mpi::environment env(argc, argv);
  mpi::communicator world;
  if(world.rank() == 0) cout << "parallel genetic-algorithm simulation" << endl;
#endif

  string confFileName;
  string dumpFileName;
  for(int i = 1; i < argc; ++i) {
    if(strcmp(argv[i], "-i") == 0) confFileName = argv[++i];
    if(strcmp(argv[i], "-f") == 0) dumpFileName = argv[++i];
  }

  ifstream confFile(confFileName.c_str());
  ifstream dumpFile(dumpFileName.c_str());

  vector<int> reorder = gaorder(confFile, dumpFile);

#ifndef SERIAL
  if(world.rank() == 0)
#endif
  {
    int n = reorder.size() - 1;
    cout << "#################### TTNS REORDER FORMAT ####################" << endl;
    for(int i = 0; i < n; ++i)
    cout << setw(3) << reorder[i] << " " << endl;
    cout << setw(3) << reorder[n] << " " << endl;
    cout << endl;

    cout << "#################### DMRG REORDER FORMAT ####################" << endl;
    for(int i = 0; i < n; ++i)
    cout << reorder[i] + 1 << ",";
    cout << reorder[n] + 1 << endl;
    cout << endl;
  }

  return 0;
}
