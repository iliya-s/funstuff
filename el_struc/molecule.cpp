#include "integrals.h"
#include "molecule.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

void molecule::build(std::string ao_integrals, std::string metric_file, std::string nuc)
{
    read_ao_integrals(ao_integrals, metric_file, nuc, ao_i1, ao_i2, metric, core_e, nuc_e, norbs, nelec, nalpha, nbeta, sz, irrep);
    built = true;
}

void molecule::rhf()
{
    if (!built)
    {
        std::cout << "Build molecule with member function molecule.build()" << std::endl;
        return;
    }

    std::cout << "Nuclear Energy: " << nuc_e << std::endl << std::endl;

    //fock matrix variables
    Eigen::MatrixXd fock = Eigen::MatrixXd::Zero(norbs, norbs);
    Eigen::MatrixXd h = fock;
    Eigen::MatrixXd v = fock;
    Eigen::MatrixXd s = fock;
    Eigen::MatrixXd d = fock;

    //build s and h matrix
    for (int i = 0; i < norbs; i++)
    {
        for (int j = 0; j < norbs; j++)
        {
            s(i,j) = metric(i,j);
            h(i,j) = ao_i1(i,j);
        }
    }
    fock = h + v;
    std::cout << "H_core: " << std::endl;
    std::cout << fock << std::endl << std::endl;;
    std::cout << "Metric: " << std::endl;
    std::cout << s << std::endl << std::endl;;

    //canonical tranformation
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);

    Eigen::MatrixXd x;
    if (es.eigenvalues().norm() == 0)
    {
        x = Eigen::MatrixXd::Identity(norbs,norbs);
    }
    else
    {
        std::vector<int> cols;
        for (int m = 0; m < norbs; m++)
        {
            if (fabs(es.eigenvalues()(m)) > 1.e-8)
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

    std::cout << "Transformation Matrix: " << std::endl;
    std::cout << x << std::endl <<std::endl;

    //inital density and orbitals
    d = Eigen::MatrixXd::Zero(norbs, norbs);
    Eigen::MatrixXd old_d = d;
    orbitals.setZero(norbs, norbs); 
    double E_total, Eold;

    std::cout << " ############### Begin Loop ################ " << std::endl << std::endl;

    for (int iter = 0; iter < 1000; iter++)
    {
        std::cout << " ############### " << iter << " ################ " << std::endl << std::endl;
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
                        v(i,j) += d(a,b) * (ao_i2(i,j,a,b) - 0.5 * ao_i2(i,b,a,j));
                    }
                }
            }
        }
        fock = h + v;
        /*
        std::cout << "Electron V Matrix: " << std::endl;
        std::cout << v << std::endl << std::endl;
        std::cout << "Fock Matrix: " << std::endl;
        std::cout << fock << std::endl << std::endl;
        */

        //diagonalize fock matrix
        Eigen::MatrixXd fock_prime = x.adjoint() * fock * x;
        /*std::cout << "Fock' Matrix: " << std::endl;
        std::cout << fock_prime << std::endl << std::endl;*/
        es.compute(fock_prime);
        orbitals_e = es.eigenvalues();
        orbitals = x * es.eigenvectors();

        //build density
        d.setZero();
        for (int i = 0; i < nalpha; i++)
        {
            d += orbitals.col(i) * orbitals.col(i).adjoint();
        }
        for (int i = 0; i < nbeta; i++)
        {
            d += orbitals.col(i) * orbitals.col(i).adjoint();
        }
        d = 0.5 * d + 0.5 * old_d;

        //calculate HF energy
        E = 0.0;
        for (int a = 0; a < norbs; a++)
        {
            for (int b = 0; b < norbs; b++)
            {
                E += 0.5 * d(a,b) * (h(a,b) + fock(a,b));
            }
        }
        E_total = E + nuc_e;

        //print data
        std::cout << "HF Electronic Energy: " << std::endl;
        std::cout << E << std::endl << std::endl;
        std::cout << "HF Total Energy: " << std::endl;
        std::cout << E_total << std::endl << std::endl;
        std::cout << "Orbital Energies: " << std::endl;
        std::cout << orbitals_e << std::endl << std::endl;;
        std::cout << "Orbitals: " << std::endl;
        std::cout << orbitals << std::endl << std::endl;;
        std::cout << "Density: " << std::endl;
        std::cout << d << std::endl << std::endl;

        //if converged
        if (fabs(E - Eold) < 1.e-5)
        {
            break;
        }
        Eold = E;
        old_d = d;
    }
}


void molecule::uhf()
{
    if (!built)
    {
        std::cout << "Build molecule with member function molecule.build()" << std::endl;
        return;
    }

    //fock matrix variables
    Eigen::MatrixXd alpha_fock = Eigen::MatrixXd::Zero(norbs, norbs);
    Eigen::MatrixXd beta_fock = alpha_fock;
    Eigen::MatrixXd h = alpha_fock;
    Eigen::MatrixXd s = alpha_fock;
    Eigen::MatrixXd total_d = alpha_fock;
    Eigen::MatrixXd alpha_d = alpha_fock;
    Eigen::MatrixXd old_alpha_d = alpha_fock;
    Eigen::MatrixXd beta_d = alpha_fock;
    Eigen::MatrixXd old_beta_d = alpha_fock;
    Eigen::MatrixXd spin_d = alpha_fock;


    //build s and h matrix
    for (int i = 0; i < norbs; i++)
    {
        for (int j = 0; j < norbs; j++)
        {
            s(i,j) = metric(i,j);
            h(i,j) = ao_i1(i,j);
        }
    }
    std::cout << "H_core: " << std::endl;
    std::cout << h << std::endl << std::endl;;
    std::cout << "Metric: " << std::endl;
    std::cout << s << std::endl << std::endl;;

    //canonical tranformation
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);

    Eigen::MatrixXd x;
    if (es.eigenvalues().norm() == 0)
    {
        x = Eigen::MatrixXd::Identity(norbs,norbs);
    }
    else
    {
        std::vector<int> cols;
        for (int m = 0; m < norbs; m++)
        {
            if (fabs(es.eigenvalues()(m)) > 1.e-8)
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

    std::cout << "Transformation Matrix: " << std::endl;
    std::cout << x << std::endl <<std::endl;

    //initalize orbitals
    //alpha_orbitals.setZero(norbs, norbs); 
    alpha_orbitals = Eigen::MatrixXd::Identity(norbs,norbs);
    beta_orbitals.setZero(norbs, norbs);
    double E_total, Eold;

    std::cout << " ############### Begin Loop ################ " << std::endl << std::endl;

    for (int iter = 0; iter < 1000; iter++)
    {
        std::cout << " ############### " << iter << " ################ " << std::endl << std::endl;
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
                        alpha_fock(i,j) += total_d(a,b) * ao_i2(i,j,a,b) - alpha_d(a,b) * ao_i2(i,b,a,j);
                        beta_fock(i,j) += total_d(a,b) * ao_i2(i,j,a,b) - beta_d(a,b) * ao_i2(i,b,a,j);
                    }
                }
            }
        }
        alpha_fock += h;
        beta_fock += h;

        //diagonalize fock matrix
        Eigen::MatrixXd alpha_fock_prime = x.adjoint() * alpha_fock * x;
        Eigen::MatrixXd beta_fock_prime = x.adjoint() * beta_fock * x;
        es.compute(alpha_fock_prime);
        alpha_orbitals_e = es.eigenvalues();
        alpha_orbitals = x * es.eigenvectors();
        es.compute(beta_fock_prime);
        beta_orbitals_e = es.eigenvalues();
        beta_orbitals = x * es.eigenvectors();

        //build density
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
        alpha_d = 0.5 * alpha_d + 0.5 * old_alpha_d;
        beta_d = 0.5 * beta_d + 0.5 * old_beta_d;
        total_d = alpha_d + beta_d;
        spin_d = alpha_d - beta_d;

        //calculate HF energy
        E = 0.0;
        for (int a = 0; a < norbs; a++)
        {
            for (int b = 0; b < norbs; b++)
            {
                E += 0.5 * (total_d(a,b) * h(a,b) + alpha_d(a,b) * alpha_fock(a,b) + beta_d(a,b) * beta_fock(a,b));
            }
        }
        E_total = E + nuc_e;

        //print data
        std::cout << "HF Electronic Energy: " << std::endl;
        std::cout << E << std::endl << std::endl;
        std::cout << "HF Total Energy: " << std::endl;
        std::cout << E_total << std::endl << std::endl;
        std::cout << "Alpha Orbital Energies: " << std::endl;
        std::cout << alpha_orbitals_e << std::endl << std::endl;;
        std::cout << "Alpha Orbitals: " << std::endl;
        std::cout << alpha_orbitals << std::endl << std::endl;;
        std::cout << "Beta Orbital Energies: " << std::endl;
        std::cout << beta_orbitals_e << std::endl << std::endl;;
        std::cout << "Beta Orbitals: " << std::endl;
        std::cout << beta_orbitals << std::endl << std::endl;;
        std::cout << "Total Density: " << std::endl;
        std::cout << total_d << std::endl << std::endl;
        std::cout << "Spin Density: " << std::endl;
        std::cout << spin_d << std::endl << std::endl;

        //if converged
        if (fabs(E - Eold) < 1.e-5)
        {
            break;
        }
        Eold = E;
        old_alpha_d = alpha_d;
        old_beta_d = beta_d;
    }
}
    
