#include <iostream>
#include "HartreeFock.h"
#include "Molecule.h"
#include "Integrals.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

void CanonicalTransform(const Eigen::MatrixXd &S, Eigen::MatrixXd &X)
{
    int dim = S.cols();
    if (dim != S.rows())
    {
        std::cout << "Input matrix not square" << std::endl;
        exit(1);
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    if (es.eigenvalues().norm() == 0)
    {
        X = Eigen::MatrixXd::Identity(dim, dim);
    }
    else
    {
        std::vector<int> cols;
        for (int m = 0; m < dim; m++)
        {
            if (std::abs(es.eigenvalues()(m)) > 1.e-10)
            {
                cols.push_back(m);
            }
        }
        X.resize(dim, cols.size());
        for (int m = 0; m < cols.size(); m++)
        {
            int index = cols[m];
            double eigval = es.eigenvalues()(index);
            X.col(m) = es.eigenvectors().col(index) / std::sqrt(eigval);
        }
    }
}

namespace HartreeFock
{
    Restricted::Restricted(Molecule &mol)
    {
        run(mol);
    }
    
    void Restricted::run(Molecule &mol)
    {
        if (!(mol.built))
        {
            std::cout << "Build molecule before running hartree fock" << std::endl;
            return;
        }
    
        //variables
        int norbs = mol.norbs;
        int nalpha = mol.nalpha;
        int nbeta = mol.nbeta;
        double core_e = mol.core_e;
        double Eold = 0.0;
        Eigen::MatrixXd fock = Eigen::MatrixXd::Zero(norbs, norbs);
        Eigen::MatrixXd density = Eigen::MatrixXd::Zero(norbs, norbs);
    
        //build s and h matrix
        Eigen::Map<Eigen::MatrixXd> h(mol.I1.store.data(), norbs, norbs);
        Eigen::Map<Eigen::MatrixXd> s(mol.S.store.data(), norbs, norbs);
    
        //canonical tranformation
        Eigen::MatrixXd x;
        CanonicalTransform(s, x);
    
        for (int iter = 0; iter < 1000; iter++)
        {
            std::cout << "____________________" << iter << "____________________" << std::endl;
            //build fock matrix
            fock = h;
            for (int i = 0; i < norbs; i++)
            {
                for (int j = 0; j < norbs; j++)
                {
                    for (int a = 0; a < norbs; a++)
                    {
                        for (int b = 0; b < norbs; b++)
                        {
                            fock(i, j) += density(a, b) * (mol.I2(i, j, a, b) - 0.5 * mol.I2(i, b, a, j));
                        }
                    }
                }
            }
            std::cout << "Fock Marix" << std::endl;
            std::cout << fock << std::endl << std::endl;
    
            //transform and diagonalize fock matrix
            Eigen::MatrixXd fock_prime = x.adjoint() * fock * x;
            std::cout << "Fock' Matrix" << std::endl;
            std::cout << fock_prime << std::endl << std::endl;
    
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(fock_prime);
            orbitals_e = es.eigenvalues();
            orbitals = x * es.eigenvectors();
            std::cout << "Orbital Energies" << std::endl;
            std::cout << orbitals_e << std::endl;
            std::cout << "Orbitals" << std::endl;
            std::cout << orbitals << std::endl << std::endl;
    
            //build density
            density.setZero(norbs, norbs);
            for (int i = 0; i < nalpha; i++)
            {
                density += orbitals.col(i) * orbitals.col(i).adjoint();
            }
            for (int i = 0; i < nbeta; i++)
            {
                density += orbitals.col(i) * orbitals.col(i).adjoint();
            }
            
            std::cout << "MO Fock Matrix" << std::endl;
            Eigen::MatrixXd MOfock = orbitals.adjoint() * fock * orbitals;
            std::cout << MOfock << std::endl << std::endl;
    
            /*
            std::cout << "MO Density Matrix" << std::endl;
            Eigen::MatrixXd MOdensity = orbitals.adjoint() * density * orbitals;
            std::cout << MOdensity << std::endl << std::endl;
            */
            std::cout << "Population Matrix" << std::endl;
            Eigen::MatrixXd N = density * s;;
            std::cout << N << std::endl << std::endl;
    
    
            //calculate HF energy
            /*
            E = 0.0;
            for (int i = 0; i < norbs; i++)
            {
                for (int j = 0; j < norbs; j++)
                {
                    E += 0.5 * density(i, j) * (h(i, j) + fock(i, j));
                }
            }
            */
            E = 0.5 * (density.array() * (h.array() + fock.array())).sum();

            std::cout << "HF Electronic Energy" << std::endl;
            std::cout << E << std::endl << std::endl;
            E += core_e;
            std::cout << "HF Energy" << std::endl;
            std::cout << E << std::endl << std::endl;
    
            if (std::abs(E - Eold) < 1.e-6)
            {
                break;
            }
            Eold = E;
        }
    } //run RHF
    
     
    Unrestricted::Unrestricted(Molecule &mol)
    {
        run(mol);
    }
    
    void Unrestricted::run(Molecule &mol)
    {
        if (!(mol.built))
        {
            std::cout << "Build molecule before running hartree fock" << std::endl;
            return;
        }
    
        //variables
        int norbs = mol.norbs;
        int nalpha = mol.nalpha;
        int nbeta = mol.nbeta;
        double core_e = mol.core_e;
        double Eold = 0.0;
        Eigen::MatrixXd spin_d = Eigen::MatrixXd::Zero(norbs, norbs); 
        Eigen::MatrixXd total_d = Eigen::MatrixXd::Zero(norbs, norbs);
        std::array<Eigen::MatrixXd, 2> fock;
        std::array<Eigen::MatrixXd, 2> d;
        d[0].setZero(norbs, norbs);
        d[1].setZero(norbs, norbs);
    
        //build s and h matrix
        Eigen::Map<Eigen::MatrixXd> h(mol.I1.store.data(), norbs, norbs);
        Eigen::Map<Eigen::MatrixXd> s(mol.S.store.data(), norbs, norbs);
    
        //canonical tranformation
        Eigen::MatrixXd x;
        CanonicalTransform(s, x);
        std::cout << x << std::endl;
    
        for (int iter = 0; iter < 1000; iter++)
        {
            std::cout << "____________________" << iter << "____________________" << std::endl;
            //build fock matrix
            fock[0] = h;
            fock[1] = h;
            for (int i = 0; i < norbs; i++)
            {
                for (int j = 0; j < norbs; j++)
                {
                    for (int a = 0; a < norbs; a++)
                    {
                        for (int b = 0; b < norbs; b++)
                        {
                            fock[0](i, j) += total_d(a, b) * mol.I2(i, j, a, b) - d[0](a, b) * mol.I2(i, b, a, j);
                            fock[1](i, j) += total_d(a, b) * mol.I2(i, j, a, b) - d[1](a, b) * mol.I2(i, b, a, j);
                        }
                    }
                }
            }
            std::cout << "Alpha Fock Marix" << std::endl;
            std::cout << fock[0] << std::endl << std::endl;
            std::cout << "Beta Fock Matrix" << std::endl;
            std::cout << fock[1] << std::endl << std::endl;
    
            //transform and diagonalize fock matrix
            std::array<Eigen::MatrixXd, 2> fock_prime;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
            for (int i = 0; i < 2; i++)
            {
                fock_prime[i] = x.adjoint() * fock[i] * x;
                es.compute(fock_prime[i]);
                orbitals_e[i] = es.eigenvalues();
                orbitals[i] = x * es.eigenvectors();
            }
            std::cout << "Alpha Fock' Matrix" << std::endl;
            std::cout << fock_prime[0] << std::endl << std::endl;
            std::cout << "Beta Fock' Matrix" << std::endl;
            std::cout << fock_prime[1] << std::endl << std::endl;
            std::cout << "Alpha Orbital Energies" << std::endl;
            std::cout << orbitals_e[0] << std::endl;
            std::cout << "Alpha Orbitals" << std::endl;
            std::cout << orbitals[0] << std::endl << std::endl;
            std::cout << "Beta Orbital Energies" << std::endl;
            std::cout << orbitals_e[1] << std::endl;
            std::cout << "Beta Orbitals" << std::endl;
            std::cout << orbitals[1] << std::endl << std::endl;
    
            /*
            std::array<Eigen::MatrixXd, 2> MOfock;
            MOfock[0] = orbitals[0].adjoint() * fock[0] * orbitals[0];
            MOfock[1] = orbitals[1].adjoint() * fock[1] * orbitals[1];
            std::cout << "Alpha MO Fock Matrix" << std::endl;
            std::cout << MOfock[0] << std::endl << std::endl;
            std::cout << "Beta MO Fock Matrix" << std::endl;
            std::cout << MOfock[1] << std::endl << std::endl;
            */
    
    
            //build density
            d[0].setZero(norbs, norbs);
            d[1].setZero(norbs, norbs);
            for (int i = 0; i < nalpha; i++)
            {
                d[0] += orbitals[0].col(i) * orbitals[0].col(i).adjoint();
            }
            for (int i = 0; i < nbeta; i++)
            {
                d[1] += orbitals[1].col(i) * orbitals[1].col(i).adjoint();
            }
            total_d = d[0] + d[1];
            spin_d = d[0] - d[1];
            std::cout << "Total Density: " << std::endl;
            std::cout << total_d << std::endl << std::endl;
            std::cout << "Spin Density: " << std::endl;
            std::cout << spin_d << std::endl << std::endl;
    
    
            //calculate HF energy
            /*
            E = 0.0;
            for (int a = 0; a < norbs; a++)
            {
                for (int b = 0; b < norbs; b++)
                {
                    E += 0.5 * (total_d(a, b) * h(a, b) + d[0](a, b) * fock[0](a, b) + d[1](a, b) * fock[1](a, b));
                }
            }
            */
            E = 0.5 * (total_d.array() * h.array() + d[0].array() * fock[0].array() + d[1].array() * fock[1].array()).sum();
            std::cout << "HF Electronic Energy" << std::endl;
            std::cout << E << std::endl << std::endl;
            E += core_e;
            std::cout << "HF Energy" << std::endl;
            std::cout << E << std::endl << std::endl;
    
            if (std::abs(E - Eold) < 1.e-6)
            {
                break;
            }
            Eold = E;
        }
    } //run UHF
}
