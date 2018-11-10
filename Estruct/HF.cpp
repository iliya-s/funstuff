#include "HF.h"
#include "Molecule.h"
#include "Integral.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;
using namespace boost;
using namespace Eigen;

void canonical_transform(MatrixXd &x, MatrixXd &s)
{
    int norbs = s.cols();
    SelfAdjointEigenSolver<MatrixXd> es(s);
    if (es.eigenvalues().norm() == 0)
    {
        x = MatrixXd::Identity(norbs, norbs);
    }
    else
    {
        vector<int> cols;
        for (int m = 0; m < norbs; m++)
        {
            if (abs(es.eigenvalues()(m)) > 1.e-8)
            {
                cols.push_back(m);
            }
        }
        x.resize(norbs, cols.size());
        for (int m = 0; m < cols.size(); m++)
        {
            int index = cols[m];
            double eigval = es.eigenvalues()(index);
            x.col(m) = es.eigenvectors().col(index) * pow(eigval, -0.5);
        }
    }
}

RHF::RHF(Molecule &mol)
{
    run(mol);
}

int RHF::run(Molecule &mol)
{
    if (!(mol.built))
    {
        cout << "Build molecule before running hartree fock" << endl;
        return 1;
    }

    //variables
    int norbs = mol.norbs;
    int nalpha = mol.nalpha;
    int nbeta = mol.nbeta;
    double nuc_e = mol.nuc_e;
    double Eold = 0.0;
    MatrixXd zero = MatrixXd::Zero(norbs,norbs);
    MatrixXd fock = zero;
    MatrixXd h = zero;
    MatrixXd v = zero;
    MatrixXd s = zero;
    vector<MatrixXd> density;
    for (int i = 0; i < 5; i++)
    {
        density.push_back(zero);
    }

    //build s and h matrix
    for (int i = 0; i < norbs; i++)
    {
        for (int j = 0; j < norbs; j++)
        {
            s(i,j) = mol.M(i,j);
            h(i,j) = mol.AI1(i,j);
        }
    }

    //canonical tranformation
    MatrixXd x = zero;
    canonical_transform(x, s);

    for (int iter = 0; iter < 1000; iter++)
    {
        cout << "____________________" << iter << "____________________" << endl;
        int d_index = iter % 5;
        int next_d_index = (iter + 1) % 5;
        //build fock matrix
        fock.setZero();
        v.setZero();
        for (int i = 0; i < norbs; i++)
        {
            for (int j = 0; j < norbs; j++)
            {
                for (int a = 0; a < norbs; a++)
                {
                    for (int b = 0; b < norbs; b++)
                    {
                        v(i,j) += density[d_index](a,b) * (mol.AI2(i,j,a,b) - 0.5 * mol.AI2(i,b,a,j));
                    }
                }
            }
        }
        fock = h + v;
        cout << "Fock Marix" << endl;
        cout << fock << endl << endl;

        //transform and diagonalize fock matrix
        MatrixXd fock_prime = x.adjoint() * fock * x;
        cout << "Fock' Matrix" << endl;
        cout << fock_prime << endl << endl;

        SelfAdjointEigenSolver<MatrixXd> es(fock_prime);
        orbitals_e= es.eigenvalues();
        orbitals = x * es.eigenvectors();
        cout << "Orbital Energies" << endl;
        cout << orbitals_e << endl;
        cout << "Orbitals" << endl;
        cout << orbitals << endl<< endl;

        //build density
        density[next_d_index].setZero();
        for (int i = 0; i < nalpha; i++)
        {
            density[next_d_index] += orbitals.col(i) * orbitals.col(i).adjoint();
        }
        for (int i = 0; i < nbeta; i++)
        {
            density[next_d_index] += orbitals.col(i) * orbitals.col(i).adjoint();
        }

        //calculate HF energy
        E = 0.0;
        for (int i = 0; i < norbs; i++)
        {
            for (int j = 0; j < norbs; j++)
            {
                E += 0.5 * density[next_d_index](i,j) * (h(i,j) + fock(i,j));
            }
        }
        cout << "HF Electronic Energy" << endl;
        cout << E << endl << endl;
        E += nuc_e;
        cout << "HF Energy" << endl;
        cout << E << endl << endl;

        if (abs(E - Eold) < 1.e-10)
        {
            break;
        }
        Eold = E;
    }
    return 0;
}



UHF::UHF(Molecule &mol)
{
    run(mol);
}

int UHF::run(Molecule &mol)
{
    if (!(mol.built))
    {
        cout << "Build molecule before running hartree fock" << endl;
        return 1;
    }

    //variables
    int norbs = mol.norbs;
    int nalpha = mol.nalpha;
    int nbeta = mol.nbeta;
    double nuc_e = mol.nuc_e;
    double Eold = 0.0;
    MatrixXd zero = MatrixXd::Zero(norbs,norbs);
    MatrixXd alpha_fock = zero;
    MatrixXd beta_fock = zero;
    MatrixXd h = zero;
    MatrixXd s = zero;
    MatrixXd alpha_d = zero;
    MatrixXd beta_d = zero;
    MatrixXd spin_d = zero;
    vector<MatrixXd> density;
    for (int i = 0; i < 5; i++)
    {
        density.push_back(zero);
    }

    //build s and h matrix
    for (int i = 0; i < norbs; i++)
    {
        for (int j = 0; j < norbs; j++)
        {
            s(i,j) = mol.M(i,j);
            h(i,j) = mol.AI1(i,j);
        }
    }

    //canonical tranformation
    MatrixXd x = zero;
    canonical_transform(x, s);

    for (int iter = 0; iter < 1000; iter++)
    {
        cout << "____________________" << iter << "____________________" << endl;
        int d_index = iter % 5;
        int next_d_index = (iter + 1) % 5;
        //build fock matrix
        alpha_fock.setZero();
        beta_fock.setZero();
        for (int i = 0; i < norbs; i++)
        {
            for (int j = 0; j < norbs; j++)
            {
                for (int a = 0; a < norbs; a++)
                {
                    for (int b = 0; b < norbs; b++)
                    {
                        alpha_fock(i,j) += density[d_index](a,b) * mol.AI2(i,j,a,b) - alpha_d(a,b) * mol.AI2(i,b,a,j);
                        beta_fock(i,j) += density[d_index](a,b) * mol.AI2(i,j,a,b) - beta_d(a,b) * mol.AI2(i,b,a,j);
                    }
                }
            }
        }
        alpha_fock += h;
        beta_fock += h;
        cout << "Alpha Fock Marix" << endl;
        cout << alpha_fock << endl << endl;
        cout << "Beta Fock Matrix" << endl;
        cout << beta_fock << endl << endl;

        //transform and diagonalize fock matrix
        MatrixXd alpha_fock_prime = x.adjoint() * alpha_fock * x;
        MatrixXd beta_fock_prime = x.adjoint() * beta_fock * x;
        cout << "Alpha Fock' Matrix" << endl;
        cout << alpha_fock_prime << endl << endl;
        cout << "Beta Fock' Matrix" << endl;
        cout << beta_fock_prime << endl << endl;

        SelfAdjointEigenSolver<MatrixXd> es(alpha_fock_prime);
        alpha_orbitals_e = es.eigenvalues();
        alpha_orbitals = x * es.eigenvectors();
        es.compute(beta_fock_prime);
        beta_orbitals_e= es.eigenvalues();
        beta_orbitals = x * es.eigenvectors();
        cout << "Alpha Orbital Energies" << endl;
        cout << alpha_orbitals_e << endl;
        cout << "Alpha Orbitals" << endl;
        cout << alpha_orbitals << endl << endl;
        cout << "Beta Orbital Energies" << endl;
        cout << beta_orbitals_e << endl;
        cout << "Beta Orbitals" << endl;
        cout << beta_orbitals << endl << endl;

        //build density
        density[next_d_index].setZero();
        alpha_d.setZero();
        beta_d.setZero();
        for (int i = 0; i < nalpha; i++)
        {
            alpha_d += alpha_orbitals.col(i) * alpha_orbitals.col(i).adjoint();
        }
        for (int i = 0; i < nbeta; i++)
        {
            beta_d += beta_orbitals.col(i) * beta_orbitals.col(i).adjoint();
        }
        density[next_d_index] = alpha_d + beta_d;
        spin_d = alpha_d - beta_d;
        cout << "Total Density: " << endl;
        cout << density[next_d_index] << endl << endl;
        cout << "Spin Density: " << endl;
        cout << spin_d << endl << endl;


        //calculate HF energy
        E = 0.0;
        for (int a = 0; a < norbs; a++)
        {
            for (int b = 0; b < norbs; b++)
            {
                E += 0.5 * (density[next_d_index](a,b) * h(a,b) + alpha_d(a,b) * alpha_fock(a,b) + beta_d(a,b) * beta_fock(a,b));
            }
        }
        cout << "HF Electronic Energy" << endl;
        cout << E << endl << endl;
        E += nuc_e;
        cout << "HF Energy" << endl;
        cout << E << endl << endl;

        if (abs(E - Eold) < 1.e-10)
        {
            break;
        }
        Eold = E;
    }
    return 0;
}
