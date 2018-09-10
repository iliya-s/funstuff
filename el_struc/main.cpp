#include "integrals.h"
#include "molecule.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

int main(int argc, char *argv[])
{
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
#endif
    std::string ao_integrals = "fcidump";
    std::string metric = "metric";
    std::string nuc = "e_nuc";

    molecule mol;
    mol.build(ao_integrals, metric, nuc);
    mol.uhf();
    //mol.rhf();
}
