#include <iostream>
#include <fstream>
#include "Determinant.h"

#include <bitset>

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
    int a = 50;
    std::cout << std::bitset<64>(a) << std::endl;
    std::cout << CountSetBits(a) << std::endl;
    std::cout << std::bitset<64>(MAXLONG) << std::endl;
    std::cout << CountSetBits(MAXLONG) << std::endl;
    cout << endl;

    long one = 1;
    long test = one << 3;
    cout << bitset<64>(test) << endl;
    cout << bitset<64>(test - one) << endl;

    Determinant c;
    cout << c(3, 0) << endl;
    c.set(3, 0, true);
    cout << c(3, 0) << endl;
    cout << c << endl;
    cout << endl;

    Determinant D, E;
    cout << D << endl;
    D.HartreeFock();
    E = D;
    cout << (E * D) << " " << (E == D) << endl;
    cout << endl;

    D.Write("hi");
    Determinant A;
    A.Read("hi");
    cout << A << endl;
    cout << endl;

    cout << A.CountSetOrbsTo(3,0) << endl;
    cout << A.Parity(3,0) << endl;
    cout << A.CountSetOrbsTo(7) << endl;
    cout << A.Parity(7) << endl;
    cout << endl;
    A *= 3.0;
    cout << A << endl << endl;
    cout << 2.0 * A << endl;
    cout << A * 2.0 << endl;
    cout << A * E << endl;
}
