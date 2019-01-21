#include "AtomicIntegrals.h"
#include "Molecule.h"
#include "HF.h"

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
    string AOintegrals = "AOfcidump";
    string metric = "metric";
    string nuc = "e_nuc";

    //Molecule mol(AOintegrals, metric, nuc);
    Molecule mol;
    mol.build(AOintegrals, metric, nuc);
    
    //RHF R(mol);
    RHF R;
    R.run(mol);

    //UHF U(mol);
}
