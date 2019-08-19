#include <iostream>
#include "Determinant.h"
#include "FundamentalOperators.h"
#include "NumberOperators.h"
#include "ExcitationOperators.h"
#include "FockVector.h"
#include "Hamiltonian.h"
#include <Eigen/Dense>

#include <bitset>
#include <map>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;

int Determinant::Nbeta;
int Determinant::Nalpha;
int Determinant::Norb;
int Determinant::Len;

int main(int argc, char *argv[])
{
    InitDetVars(0, 10, 8);
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
#endif 
    {
        long r = 2;
        cout << r << " in bits: " << bitset<64>(r) << endl;
        cout << -r << " in bits: " << bitset<64>(-r) << endl;
        long f = 194982;
        cout << f << " in bits: " << bitset<64>(f) << endl;
        cout << -f << " in bits: " << bitset<64>(-f) << endl;
        cout << f << " xor " << -f << ": " << bitset<64>(f & -f) << endl;
        cout << f << " xor (" << f << " and " << -f << "):" << bitset<64>(f ^ (f & -f)) << endl;
        cout << "__builtin_ctzl(" << f << ") = " << __builtin_ctzl(f) << endl;
    }

    {
        cout << "Determinant" << endl;
        cout << "set orbital 3alpha" << endl;
        Determinant c;
        cout << c(3, 0) << endl;
        c.set(3, 0, true);
        cout << c(3, 0) << endl;
        cout << c << endl;
        cout << endl;

        cout << "Hartree Fock" << endl;
        Determinant D, M;
        cout << D << endl;
        D.HartreeFock();
        cout << D << endl;
        M = D;
        cout << "Overlap" << endl;
        cout << (M * D) << " " << (M == D) << endl;
        cout << endl;

        cout << "save and write" << endl;
        D.write("hi");
        Determinant A;
        A.read("hi");
        cout << A << endl;
        cout << endl;

        cout << "Count set orbs and parity" << endl;
        cout << A.CountSetOrbsTo(3,0) << endl;
        cout << A.parity(3,0) << endl;
        cout << A.CountSetOrbsTo(7) << endl;
        cout << A.parity(7) << endl;
        cout << endl;
        A *= 3.0;
        cout << A << endl << endl;
        cout << A / 2 << endl;
        cout << A * 2.0 << endl;
        cout << A * M << endl;
        cout << A * A << endl;
        cout << endl << endl;
    }

    {
        Determinant D;
        D.HartreeFock();
        cout << "LadderOperators" << endl;
        cout << D << endl;
        Operator::Annihilation a;
        Operator::Creation a_dag;
        cout << a(1) * D << endl;
        cout << a_dag(14) * (a(1) * D) << endl;

        cout << "NumberOperators" << endl;
        Operator::OccupationNumber n;
        Operator::ParticleNumber N;
        cout << n(2) * D << endl;
        cout << N * D << endl;

        cout << "ExcitationOperators" << endl;
        Operator::Excitation E;
        cout << "Single Excitation" << endl;
        D.HartreeFock();
        Determinant L = E(11, 8) * D;
        cout << "D " <<  D << endl;
        cout << "L " <<  L << endl;
        int i, j, o, p;
        OneDiffOrbIndices(D, L, i, j);
        cout << i << " " << j << endl;
        cout << "Double Excitation" << endl;
        L = E(15, 6) * L;
        cout << "D " <<  D << endl;
        cout << "L " <<  L << endl;
        TwoDiffOrbIndices(D, L, i, j, o, p);
        cout << i << " " << j << " " << o << " " << p << endl;
    }

    {
        std::cout << std::endl << "FCI-vector" << std::endl;
    
        std::vector<std::vector<int>> alpha, beta;
        GenerateCombinations(4, 2, alpha);
        GenerateCombinations(4, 3, beta);
        std::cout << std::endl << "alpha comb" << std::endl;
        std::cout << alpha.size() << endl << endl;
        for (int i = 0; i < alpha.size(); i++)
        {
            for (int j = 0; j < alpha[i].size(); j++)
            {
                std::cout << alpha[i][j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << "beta comb" << std::endl;
        std::cout << beta.size() << endl << endl;
        for (int i = 0; i < beta.size(); i++)
        {
            for (int j = 0; j < beta[i].size(); j++)
            {
                std::cout << beta[i][j];
            }
            std::cout << std::endl;
        }
        cout << endl;
    
        FCIVector V;
        cout << V.size() << endl;
    }
    
    {
        cout << "Hamiltonian" << std::endl;
        cout << "hf det and energy" << std::endl;
        Operator::Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        FCIVector FCI;
        cout << "size of fock space: " << FCI.size() << endl;
        cout << "Diagonal" << std::endl;
        Eigen::VectorXd V;
        H.diagonal(FCI, V);
        cout << V.transpose() << endl << endl;
        cout << "Matrix" << std::endl;
        Eigen::MatrixXd Ham;
        H.matrix(FCI, Ham);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ham);
        if (Ham.cols() <= 10)
        {
            cout << Ham << endl << endl;
            cout << "\nEigenproblem Solution" << std::endl;
            std::cout << es.eigenvalues().transpose() << std::endl << std::endl;
            std::cout << es.eigenvectors() << std::endl;
        }
        cout << "\nGround state" << std::endl;
        cout << "Energy: " << es.eigenvalues()(0) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        FCI.update(gs);
        Eigen::VectorXd test;
        FCI.vector(test);
        cout << test.transpose() * gs << endl;
        FCI.trim();
        std::cout << FCI << std::endl;
    }
}
