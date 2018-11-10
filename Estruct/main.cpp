#include "Integral.h"
#include "Molecule.h"
#include "HF.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
using namespace boost;
#endif

using namespace std;

int main(int argc, char *argv[])
{
#ifndef SERIAL
    mpi::environment env(argc, argv);
    mpi::communicator world;
#endif
    string ao_integrals = "fcidump";
    string metric = "metric";
    string nuc = "e_nuc";

    Molecule mol(ao_integrals, metric, nuc);
    //mol.build(ao_integrals, metric, nuc);
    
    /*
    RHF r;
    r.run(mol);
    */

    UHF u(mol);

    /*
    mol.integral_transform("rhf");
    double E = 0.0;

    cout << mol.mo_i1.store << endl;

    cout << mol.mo_i2(0,0,0,0) << endl;
    cout << mol.mo_i2(1,1,1,1) << endl;
    cout << mol.mo_i2(1,1,0,0) << endl;
    cout << mol.mo_i2(0,1,0,1) << endl;
    for (int i = 0; i < mol.nelec; i++)
    {
        E += mol.mo_i1(i,i);
        for (int j = 0; j < mol.nelec; j++)
        {
            E += mol.mo_i2(i,i,j,j) - mol.mo_i2(i,j,i,j);
        }
    }
    cout << E << endl;
    */
}
