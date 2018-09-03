#include "ao_integrals.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
#ifndef SERIAL
#include </boost/interprocess/shared_memory_object.hpp>
#include <boost/mpi.hpp>
#endif

int main(int argc, char *argv[])
{
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
#endif
    std::string fcidump = "fcidump";
    std::string metric = "metric";
    std::string nuc = "e_nuc";

    //integral variables
    one_int I1, S;
    two_int I2;
    double core_e, nuc_e;
    int norbs, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;

    read_ao_integrals(fcidump, metric, nuc, I1, I2, S, core_e, nuc_e, norbs, nelec, nalpha, nbeta, sz, irrep);
    std::cout << "Nuclear energy: " << nuc_e << std::endl << std::endl;

    /*
    std::cout << "2-integrals: " << std::endl;
    std::cout << I2(0,0,0,0) << " " << I2(1,1,1,1) << std::endl;
    std::cout << I2(1,0,0,0) << " " << I2(0,1,0,0) << std::endl;
    std::cout << I2(0,0,1,0) << " " << I2(0,0,0,1) << std::endl;
    std::cout << I2(1,1,0,0) << " " << I2(0,0,1,1) << std::endl;
    std::cout << I2(1,0,1,0) << " " << I2(0,1,0,1) << std::endl;
    std::cout << I2(0,1,1,0) << " " << I2(1,0,0,1) << std::endl;
    std::cout << I2(1,1,1,0) << " " << I2(1,1,0,1) << std::endl;
    std::cout << I2(1,0,1,1) << " " << I2(0,1,1,1) << std::endl << std::endl;
    */

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
            s(i,j) = S(i,j);
            h(i,j) = I1(i,j);
        }
    }
    fock = h + v;
    std::cout << "H_core: " << std::endl;
    std::cout << fock << std::endl << std::endl;;
    std::cout << "metric: " << std::endl;
    std::cout << s << std::endl << std::endl;;
    
    //canonical tranformation
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
    std::cout << "Eig for metric: " << std::endl;
    std::cout << es.eigenvalues() << std::endl;
    std::cout << es.eigenvectors() << std::endl << std::endl;
    std::vector<int> cols;
    for (int m = 0; m < norbs; m++)
    {
        if (fabs(es.eigenvalues()(m)) > 1.e-8)
        {
            cols.push_back(m);
        }
    }
    Eigen::MatrixXd x(norbs, cols.size());
    for (int m = 0; m < cols.size(); m++)
    {
        int index = cols[m];
        double eigval = es.eigenvalues()(index);
        x.col(m) = es.eigenvectors().col(index) * pow(eigval, -0.5);
    }
    /* 
    //symmetric transformation
    Eigen::MatrixXd x(norbs, norbs);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(s);
    Eigen::VectorXd s_half_inv(norbs);
    for (int i = 0; i < norbs; i++)
    {
        s_half_inv(i) = 1.0 / pow(es.eigenvalues()(i), 0.5);
    }
    x = es.eigenvectors() * s_half_inv.asDiagonal() * es.eigenvectors().adjoint();
    */
    std::cout << "Transformation matrix: " << std::endl;
    std::cout << x << std::endl <<std::endl;

    //inital density and orbitals
    //d = Eigen::MatrixXd::Identity(norbs, norbs);
    d = Eigen::MatrixXd::Zero(norbs, norbs);
    Eigen::MatrixXd orbitals(norbs,norbs);
    
    std::cout << " ###############begin loop################ " << std::endl << std::endl;
    //while(true)
    for (int p = 0; p < 6; p++)
    {

        std::cout << " ############### " << p << " ################ " << std::endl << std::endl;
        //build fock matrix
        fock.setZero();
        v.setZero();
        for (int i = 0; i < norbs; i++)
        {
            for (int j = 0; j < norbs; j++)
            {
                for (int a = 0; a < nelec; a++)
                {
                    for (int b = 0; b < nelec; b++)
                    {
                        v(i,j) += d(a,b) * (I2(i,j,a,b) - 0.5 * I2(i,b,a,j));
                    }
                }
            }
        }
        std::cout << "electron v matrix: " << std::endl;
        std::cout << v << std::endl << std::endl;
        fock = h + v;
        std::cout << "fock matrix: " << std::endl;
        std::cout << fock << std::endl << std::endl;
    
        //diagonalize fock matrix
        Eigen::MatrixXd fock_prime = x.adjoint() * fock * x;
        std::cout << "fock_prime matrix: " << std::endl;
        std::cout << fock_prime << std::endl << std::endl;
        es.compute(fock_prime);
        Eigen::VectorXd e = es.eigenvalues();
        orbitals = x * es.eigenvectors();

        d.setZero();
        for (int i = 0; i < nalpha; i++)
        {
            d += orbitals.col(i) * orbitals.col(i).adjoint();
        }
        for (int i = 0; i < nbeta; i++)
        {
            d += orbitals.col(i) * orbitals.col(i).adjoint();
        }

        double E = 0.0;
        for (int a = 0; a < norbs; a++)
        {
            for (int b = 0; b < norbs; b++)
            {
                E += 0.5 * d(a,b) * (h(a,b) + fock(a,b));
            }
        }
        double E_total = E + nuc_e;

        std::cout << "HF Electronic Energy: " << std::endl;
        std::cout << E << std::endl << std::endl;
        std::cout << "HF Total Energy: " << std::endl;
        std::cout << E_total << std::endl << std::endl;
        std::cout << "Orbital energies: " << std::endl;
        std::cout << e << std::endl << std::endl;;
        std::cout << "Orbitals: " << std::endl;
        std::cout << orbitals << std::endl << std::endl;;
        std::cout << "Density: " << std::endl;
        std::cout << d << std::endl << std::endl;
    }
}


