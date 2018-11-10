#include "Integral.h"
#include "Molecule.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

Molecule::Molecule(string ao_fcidump, string metric_file, string nuc_file)
{
    build(ao_fcidump, metric_file, nuc_file);
}

int Molecule::build(string ao_fcidump, string metric_file, string nuc_file)
{
    built = true;
    return read_ao_integrals(ao_fcidump, metric_file, nuc_file, AI1, AI2, M, core_e, nuc_e, norbs, nelec, nalpha, nbeta, sz, irrep);
}

/* 
void molecule::integral_transform(std::string hf)
{
    if (hf == "rhf")
    {
        mo_i1.store.setZero(norbs,norbs);
        mo_i1.store = orbitals.transpose() * ao_i1.store * orbitals;

        mo_i2.store.resize((norbs * norbs * norbs * norbs), 0.0);
        double *temp;
        temp = new double[norbs][norbs][norbs][norbs];
        double *temp1;
        temp1 = new double[norbs][norbs][norbs][norbs];
        double *temp2;
        temp2 = new double[norbs][norbs][norbs][norbs];
        double *temp3;
        temp3 = new double[norbs][norbs][norbs][norbs];
        
        for (int i = 0; i < norbs; i++)
        {
            for (int j = 0; j < norbs; j++)
            {
                for (int a = 0; a < norbs; a++)
                {
                    for (int b = 0; b < norbs; b++)
                    {
                        temp[i][j][a][b] = 0.0;
                        temp1[i][j][a][b] = 0.0;
                        temp2[i][j][a][b] = 0.0;
                        temp3[i][j][a][b] = 0.0;
                    }
                }
            }
        }


        for (int i = 0; i < norbs; i++)
        {
            for (int j = 0; j < norbs; j++)
            {
                for (int a = 0; a < norbs; a++)
                {
                    for (int b = 0; b < norbs; b++)
                    {
                        mo_i2(i,j,a,b) = temp[i][j][a][b];
                    }
                }
            }
        }

    }
}
*/

