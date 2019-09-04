#include <iostream>
#include <algorithm>
#include "Determinant.h"
#include "FundamentalOperators.h"
#include "NumberOperators.h"
#include "ExcitationOperators.h"
#include "FockVector.h"
#include "Hamiltonian.h"
#include "Davidson.h"
#include <Eigen/Dense>

#include <bitset>
#include <map>
//#include <unordered_map>
#include <ctime>

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
        /*
        std::vector<int> r(100000000);
        std::iota(r.begin(), r.end(), 0);
        std::random_shuffle(r.begin(), r.end());

        auto begin = std::time(nullptr);
        std::unordered_map<int, int> A;
        for (int i = 0; i < r.size(); i++)
        {
            A.insert(std::make_pair(r[i], i));
        }
        auto end = std::time(nullptr);
        cout << "insert into hashtable: " << end - begin << endl;

        begin = std::time(nullptr);
        std::vector<int> B;
        for (int i = 0; i < r.size(); i++)
        {
            B.push_back(r[i]);
        }
        end = std::time(nullptr);
        cout << "insert into vector: " << end - begin << endl;

        begin = std::time(nullptr);
        std::sort(B.begin(), B.end());
        end = std::time(nullptr);
        cout << "time to sort vector: " << end - begin << endl;

        begin = std::time(nullptr);
        for (int i = 0; i < r.size(); i++)
        {
            int val = A.at(i);
        }
        end = std::time(nullptr);
        cout << "access all elements in hashtable: " << end - begin << endl;

        begin = std::time(nullptr);
        for (int i = 0; i < r.size(); i++)
        {
            auto it = std::lower_bound(B.begin(), B.end(), i);
        }
        end = std::time(nullptr);
        cout << "access all elements in vector: " << end - begin << endl;
        */
    }
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
        auto ij = OneDiffOrbIndices(D, L);
        i = ij.first;
        j = ij.second;
        cout << i << " " << j << endl;
        cout << "Double Excitation" << endl;
        L = E(15, 6) * L;
        cout << "D " <<  D << endl;
        cout << "L " <<  L << endl;
        auto ijab = TwoDiffOrbIndices(D, L);
        i = ijab.first.first;
        j = ijab.first.second;
        o = ijab.second.first;
        p = ijab.second.second;
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
    
    /*
    {
        cout << "Hamiltonian" << std::endl;
        cout << "hf det and energy" << std::endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        FCIVector FCI;
        cout << "size of fock space: " << FCI.size() << endl;
        //cout << "Diagonal" << std::endl;
        Eigen::VectorXd V;
        H.diagonal(FCI, V);
        //cout << V.transpose() << endl << endl;
        //cout << "Matrix" << std::endl;
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
        double val = es.eigenvalues()(0);
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        Eigen::VectorXd r = Ham * gs - val * gs;
        cout << "residual: " << r.norm() << endl;
        FCI.update(gs);
        Eigen::VectorXd test;
        FCI.vector(test);
        cout << test.transpose() * gs << endl;
        FCI.trim(1.e-3);
        std::cout << FCI << std::endl;
    }
    */

    {
        cout << "generating connected determinants" << endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        std::vector<Determinant> V;
        hf.excitations(V);
        cout << V.size() << endl;
        /*
        for (int i = 0; i < V.size(); i++)
        {
            cout << V[i] << endl;
        }
        */
    }

    {
        /*
        cout << "\n\nCISD" << endl;
        cout << "hf det and energy" << std::endl;

        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        CISDVector CI;
        H.space(CI);
        cout << "size of fock space: " << CI.size() << endl;
        //cout << "Diagonal" << std::endl;
        Eigen::VectorXd V;
        H.diagonal(V);
        //cout << V.transpose() << endl << endl;
        cout << "Matrix" << std::endl;
        Eigen::MatrixXd Ham;
        H.matrix(Ham);
        auto begin = std::time(nullptr);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ham);
        auto end = std::time(nullptr);
        if (Ham.cols() <= 10)
        {
            cout << Ham << endl << endl;
            cout << "\nEigenproblem Solution" << std::endl;
            std::cout << es.eigenvalues().transpose() << std::endl << std::endl;
            std::cout << es.eigenvectors() << std::endl;
        }
        cout << "\n" << std::endl;
        double val = es.eigenvalues()(0);
        cout << "Ground state Energy: " << es.eigenvalues()(0) << endl;
        cout << "1st excited state Energy: " << es.eigenvalues()(1) << endl;
        cout << "2st excited state Energy: " << es.eigenvalues()(2) << endl;
        cout << "3rd excited state Energy: " << es.eigenvalues()(3) << endl;
        cout << "4th excited state Energy: " << es.eigenvalues()(4) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        Eigen::VectorXd r = Ham * gs - val * gs;
        cout << "residual: " << r.norm() << endl;
        CI.update(gs);
        Eigen::VectorXd test;
        CI.vector(test);
        cout << test.transpose() * gs << endl;
        CI.trim(1.e-3);
        std::cout << CI << std::endl;
        std::cout << "total time: " << end - begin << std::endl;
        */
    }

    /*
    {
        cout << "\n\nCISD with Davidson" << endl;
        cout << "hf det and energy" << std::endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        CISDVector CI;
        cout << "size of fock space: " << CI.size() << endl;
        //cout << "Diagonal" << std::endl;
        Eigen::VectorXd V;
        H.diagonal(CI, V);
        //cout << V.transpose() << endl << endl;
        cout << "Matrix" << std::endl;
        Eigen::MatrixXd Ham;
        H.matrix(CI, Ham);
        MatrixMult Mat(H, CI);
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ham);
        //Ham.diagonal().array() += 0.1;
        Davidson es;
        auto begin = std::time(nullptr);
        int numIter = es.run(Ham);
        auto end = std::time(nullptr);
        if (Ham.cols() <= 10)
        {
            cout << Ham << endl << endl;
            cout << "\nEigenproblem Solution" << std::endl;
            std::cout << es.eigenvalues().transpose() << std::endl << std::endl;
            std::cout << es.eigenvectors() << std::endl;
        }
        cout << "number of iterations: " << numIter << std::endl;
        double val = es.eigenvalues()(0);
        cout << "\nGround state" << std::endl;
        cout << "Energy: " << es.eigenvalues()(0) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        Eigen::VectorXd r = Ham * gs - val * gs;
        cout << "residual: " << r.norm() << endl;
        CI.update(gs);
        Eigen::VectorXd test;
        CI.vector(test);
        cout << test.transpose() * gs << endl;
        CI.trim(1.e-3);
        std::cout << CI << std::endl;
        std::cout << "total time: " << end - begin << std::endl;
    }
    */
    {
        cout << "\n\nCIS with direct Davidson" << endl;
        cout << "hf det and energy" << std::endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        CISVector CI;
        H.space(CI);
        cout << "size of fock space: " << CI.size() << endl;
        Davidson es;
        //DirectMatrixMult Mat(H);
        auto begin = std::time(nullptr);
        int numIter = es.run(H);
        auto end = std::time(nullptr);
        cout << "number of iterations: " << numIter << std::endl;
        double val = es.eigenvalues()(0);
        cout << "\n" << std::endl;
        cout << "Ground state Energy: " << es.eigenvalues()(0) << endl;
        //cout << "1st excited state Energy: " << es.eigenvalues()(1) << endl;
        //cout << "2st excited state Energy: " << es.eigenvalues()(2) << endl;
        //cout << "3rd excited state Energy: " << es.eigenvalues()(3) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        CI.update(gs);
        Eigen::VectorXd Hgs = H * gs;
        Eigen::VectorXd r = Hgs - val * gs;
        cout << "residual: " << r.norm() << endl;
        std::cout << "Compactness of Hamiltonian: " << ((double)H.size() / (double)(CI.size() * CI.size())) << endl;
        /*
        std::vector<Determinant> dets;
        hf.connected(dets);
        for (int i = 0; i < 10; i++)
        {
            cout << CI(dets[i]) << endl;
        }
        */
        //CI.trim(1.e-3);
        //std::cout << CI << std::endl;
        std::cout << "total time: " << end - begin << std::endl;
        cout << "\nCorrelation energy\n";
        cout << val - H.energy(hf) << endl;
        cout << "Davdison correction\n";
        cout << (1.0 - CI(hf) * CI(hf)) * (val - H.energy(hf));
    }
    {
        cout << "\n\nCID with direct Davidson" << endl;
        cout << "hf det and energy" << std::endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        CIDVector CI;
        H.space(CI);
        cout << "size of fock space: " << CI.size() << endl;
        Davidson es;
        //DirectMatrixMult Mat(H);
        auto begin = std::time(nullptr);
        int numIter = es.run(H);
        auto end = std::time(nullptr);
        cout << "number of iterations: " << numIter << std::endl;
        double val = es.eigenvalues()(0);
        cout << "\n" << std::endl;
        cout << "Ground state Energy: " << es.eigenvalues()(0) << endl;
        //cout << "1st excited state Energy: " << es.eigenvalues()(1) << endl;
        //cout << "2st excited state Energy: " << es.eigenvalues()(2) << endl;
        //cout << "3rd excited state Energy: " << es.eigenvalues()(3) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        CI.update(gs);
        Eigen::VectorXd Hgs = H * gs;
        Eigen::VectorXd r = Hgs - val * gs;
        cout << "residual: " << r.norm() << endl;
        std::cout << "Compactness of Hamiltonian: " << ((double)H.size() / (double)(CI.size() * CI.size())) << endl;
        /*
        std::vector<Determinant> dets;
        hf.connected(dets);
        for (int i = 0; i < 10; i++)
        {
            cout << CI(dets[i]) << endl;
        }
        */
        CI.trim(1.e-3);
        //std::cout << CI << std::endl;
        std::cout << "total time: " << end - begin << std::endl;
        cout << "\nCorrelation energy\n";
        cout << val - H.energy(hf) << endl;
        cout << "Davdison correction\n";
        cout << (1.0 - CI(hf) * CI(hf)) * (val - H.energy(hf));
    }

    {
        cout << "\n\nCISD with direct Davidson" << endl;
        cout << "hf det and energy" << std::endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        CISDVector CI;
        H.space(CI);
        cout << "size of fock space: " << CI.size() << endl;
        Davidson es;
        //DirectMatrixMult Mat(H);
        auto begin = std::time(nullptr);
        int numIter = es.run(H);
        auto end = std::time(nullptr);
        cout << "number of iterations: " << numIter << std::endl;
        double val = es.eigenvalues()(0);
        cout << "\n" << std::endl;
        cout << "Ground state Energy: " << es.eigenvalues()(0) << endl;
        //cout << "1st excited state Energy: " << es.eigenvalues()(1) << endl;
        //cout << "2st excited state Energy: " << es.eigenvalues()(2) << endl;
        //cout << "3rd excited state Energy: " << es.eigenvalues()(3) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        CI.update(gs);
        Eigen::VectorXd Hgs = H * gs;
        Eigen::VectorXd r = Hgs - val * gs;
        cout << "residual: " << r.norm() << endl;
        std::cout << "Compactness of Hamiltonian: " << ((double)H.size() / (double)(CI.size() * CI.size())) << endl;
        /*
        std::vector<Determinant> dets;
        hf.connected(dets);
        for (int i = 0; i < 10; i++)
        {
            cout << CI(dets[i]) << endl;
        }
        */
        CI.trim(1.e-3);
        //std::cout << CI << std::endl;
        std::cout << "total time: " << end - begin << std::endl;
        cout << "\nCorrelation energy\n";
        cout << val - H.energy(hf) << endl;
        cout << "Davdison correction\n";
        cout << (1.0 - CI(hf) * CI(hf)) * (val - H.energy(hf));
    }
    {
        cout << "\n\nFCI with direct Davidson" << endl;
        cout << "hf det and energy" << std::endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;
        FCIVector CI;
        H.space(CI);
        cout << "size of fock space: " << CI.size() << endl;
        //cout << "Diagonal" << std::endl;
        //Eigen::VectorXd V = Eigen::VectorXd::Random(CI.size());
        //H.diagonal(CI, V);
        //cout << V.transpose() << endl << endl;
        /*
        Eigen::MatrixXd Ham;
        H.matrix(CI, Ham);
        Eigen::VectorXd mult = Ham * V;
        cout << "Matrix mult" << std::endl;
        cout << mult.transpose() << endl;
        cout << "direct mult" << std::endl;
        Eigen::VectorXd direct = H.multiply(CI, V);
        cout << direct.transpose() << endl;
        */
        //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ham);
        //Ham.diagonal().array() += 0.1;
        Davidson es;
        //DirectMatrixMult Mat(H);
        auto begin = std::time(nullptr);
        int numIter = es.run(H);
        auto end = std::time(nullptr);
        cout << "number of iterations: " << numIter << std::endl;
        double val = es.eigenvalues()(0);
        cout << "\n" << std::endl;
        cout << "Ground state Energy: " << es.eigenvalues()(0) << endl;
        //cout << "1st excited state Energy: " << es.eigenvalues()(1) << endl;
        //cout << "2st excited state Energy: " << es.eigenvalues()(2) << endl;
        //cout << "3rd excited state Energy: " << es.eigenvalues()(3) << endl;
        Eigen::VectorXd gs = es.eigenvectors().col(0);
        Eigen::VectorXd Hgs = H * gs;
        Eigen::VectorXd r = Hgs - val * gs;
        cout << "residual: " << r.norm() << endl;
        CI.update(gs);
        /*
        std::vector<Determinant> dets;
        hf.connected(dets);
        for (int i = 0; i < 10; i++)
        {
            cout << CI(dets[i]) << endl;
        }
        */
        std::cout << "Compactness of Hamiltonian: " << ((double)H.size() / (double)(CI.size() * CI.size())) << endl;
        CI.trim(1.e-3);
        //std::cout << CI << std::endl;
        std::cout << "total time: " << end - begin << std::endl;
        cout << "\nCorrelation energy\n";
        cout << val - H.energy(hf) << endl;
    }
    {
        /*
        cout << "\n\nMult" << endl;
        Hamiltonian H;
        Determinant hf;
        hf.HartreeFock();
        cout << H.energy(hf) * hf << endl;

        CISDVector CI;
        H.space(CI);

        Eigen::MatrixXd Ham;
        H.matrix(Ham);

        Eigen::VectorXd V = Eigen::VectorXd::Random(CI.size());
        Eigen::VectorXd direct;
        H.Fmultiply(V, direct);

        Eigen::VectorXd mult = Ham * V;
        cout << "Matrix mult" << std::endl;
        cout << mult.transpose() << endl;
        cout << "Direct mult" << std::endl;
        cout << direct.transpose() << endl;
        */
    }
    {
        cout << endl << endl;
        Hamiltonian H;
        //Integral::OneElectron I1;
        //Integral::TwoElectron I2;
        //double core_e;
        //int norb, nelec, nalpha, nbeta, sz;
        //std::vector<int> irrep;
        //ReadFCIDUMP("FCIDUMP", I1, I2, core_e, norb, nelec, nalpha, nbeta, sz, irrep);
        //InitDetVars(sz, nelec, norb);
        //Integral::HeatBath::OneElectron HBI1(I1, I2);
        //Integral::HeatBath::TwoElectron HBI2(I1, I2);

        /*
        for (auto it = HBI1.Store.begin(); it != HBI1.Store.end(); it++)
        {
            int i = it->first;
            auto temp = it->second;
            for (auto it1 = temp.begin(); it1 != temp.end(); it1++)
            {
                float val = it1->first;
                int j = it1->second;
                cout << i << " -> " << j << " = " << val << endl;
            }
            cout << endl;
        }
        cout << endl;

        for (auto it = HBI2.Store.begin(); it != HBI2.Store.end(); it++)
        {
            auto ij = it->first;
            auto temp = it->second;
            for (auto it1 = temp.begin(); it1 != temp.end(); it1++)
            {
                float val = it1->first;
                auto ab = it1->second;
                cout << ij.first << " " << ij.second << " -> " << ab.first << " " << ab.second << " = " << val << endl;
            }
            cout << endl;
        }
        cout << endl;
        */

        Determinant hf;
        hf.HartreeFock();
        std::vector<Determinant> dets, Sdets;
        hf.excitations(dets);
        hf.screenedExcitations(H.heatBathOneElectron(), H.heatBathTwoElectron(), Sdets, 0.0);
        cout << dets.size() << " " << Sdets.size() << endl;



    }
}
