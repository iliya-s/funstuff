#include <iostream>
#include "Determinant.h"
#include "FundamentalOperators.h"
#include "NumberOperators.h"
#include "ExcitationOperators.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;

int Determinant::nbeta;
int Determinant::nalpha;
int Determinant::norbs;
int Determinant::len;

int main(int argc, char **argv)
{
    InitDetVars(0, 10, 10);
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
#endif

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
    cout << 2.0 * A << endl;
    cout << A * 2.0 << endl;
    cout << A * M << endl;
    cout << endl << endl;

    cout << "LadderOperators" << endl;
    cout << D << endl;
    Operator::Annihilation a;
    Operator::Creation a_dag;
    cout << a(1) * D << endl;
    cout << a_dag(14) * D << endl;

    cout << "NumberOperators" << endl;
    Operator::OccupationNumber n;
    Operator::ParticleNumber N;
    cout << n(2) * D << endl;
    cout << N * D << endl;

    cout << "ExcitationOperators" << endl;
    Operator::Excitation E(10, 8);
    D.HartreeFock();
    cout << D << endl;
    Determinant L = E * D;
    cout << L << endl;
    cout << E.adjoint() * L << endl;
}
