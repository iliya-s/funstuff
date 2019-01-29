#include "Molecule.h"
#include "HartreeFock.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;

int main(int argc, char **argv)
{
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
#endif
    string AOFCIDUMP = "AOFCIDUMP";
    string METRIC = "METRIC";

    Molecule mol;
    mol.build(AOFCIDUMP, METRIC);
    
    HartreeFock::Restricted R;
    R.run(mol);
}
