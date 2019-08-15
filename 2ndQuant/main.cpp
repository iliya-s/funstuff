#include <iostream>
#include "Determinant.h"
#include "FundamentalOperators.h"
#include "NumberOperators.h"
#include "ExcitationOperators.h"
#include "Hamiltonian.h"

#include <bitset>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;

int Determinant::nbeta;
int Determinant::nalpha;
int Determinant::norbs;
int Determinant::len;



//int main(int argc, char **argv)
int main(void)
{
    InitDetVars(0, 10, 10);
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
#endif
    long r = 2;
    cout << bitset<64>(r) << endl;
    cout << bitset<64>(-r) << endl;
    long f = 194982;
    cout << bitset<64>(f) << endl;
    cout << bitset<64>(-f) << endl;
    cout << bitset<64>(f & -f) << endl;
    cout << bitset<64>(f ^ (f & -f)) << endl;
    cout << __builtin_ctzl(f) << endl;

    cout << "Determinant" << endl;
    Determinant c;
    cout << c(3, 0) << endl;
    c.set(3, 0, true);
    cout << c(3, 0) << endl;
    cout << c << endl;
    cout << endl;

    Determinant D, M;
    cout << D << endl;
    D.HartreeFock();
    M = D;
    cout << (M * D) << " " << (M == D) << endl;
    cout << endl;

    D.write("hi");
    Determinant A;
    A.read("hi");
    cout << A << endl;
    cout << endl;

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
