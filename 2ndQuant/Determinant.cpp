#include "Determinant.h"
#include <fstream>

//initialize static variables
void InitDetVars(int spin, int nelec, int norbs)
{
    Determinant::norbs = norbs;
    Determinant::len = norbs / 64 + 1;
    Determinant::nalpha = nelec / 2 + spin;
    Determinant::nbeta = nelec - Determinant::nalpha;
}

//constructor
Determinant::Determinant()
{
    String = new long *[2];
    String[0] = new long [len];
    String[1] = new long [len];
    long zero = 0;
    for (int i = 0; i < len; i++)
    {
        String[0][i] = zero;
        String[1][i] = zero;
    }
}

Determinant::Determinant(const Determinant &D)
{
    for (int i = 0; i < len; i++)
    {
        String[0][i] = D.String[0][i];
        String[1][i] = D.String[1][i];
    }
}

//deconstructor
Determinant::~Determinant()
{
    delete[] String[0];
    delete[] String[1];
    delete[] String;
}

//operators
//  copy/assignment operator
Determinant &Determinant::operator=(const Determinant &RHS)
{
    for (int i = 0; i < len; i++)
    {
        String[0][i] = RHS.String[0][i];
        String[1][i] = RHS.String[1][i];
    }
    return *this;
}

//  multiplication operator overloaded for overlap = <bra|ket>
int Determinant::operator*(const Determinant &RHS) const
{
    for (int i = 0; i < len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i])
            return 0;
    }
    return 1;
}

//  equivalence comparison
bool Determinant::operator==(const Determinant &RHS) const
{
    for (int i = 0; i < len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i])
            return false;
    }
    return true;
}

//  output stream
std::ostream &operator<<(std::ostream &os, const Determinant &D)
{
    for (int i = 0, norbs = Determinant::norbs; i < norbs; i++)
    {
        bool alpha = D(i, 0);
        bool beta = D(i, 1);
        if (!alpha && !beta)
            os << "0 ";
        else if (alpha && !beta)
            os << "a ";
        else if (!alpha && beta)
            os << "b ";
        else if (alpha && beta)
            os << "2 ";
        if ((i + 1) % 4 == 0)
            os << "  ";
    }
    return os;
}

//getter
inline bool Determinant::operator()(int orbital, int spin) const
{
    long index = orbital / 64, bit = orbital % 64, one = 1;
    long test = String[spin][index] & (one << bit);
    if (test == 0)
        return false;
    else
        return true;
}

//setter
inline void Determinant::set(int orbital, int spin, bool occupancy)
{
    long index = orbital / 64, bit = orbital % 64, one = 1;
    if (occupancy)
        String[spin][index] |= (one << bit);
    else
        String[spin][index] &= ~(one << bit);
}

//Counts occupied orbitals up to but not including specified orbital
int Determinant::CountSetOrbsTo(int orbital, int spin) const
{
    long index = orbital / 64, bit = orbital % 64, one = 1;
    int count = 0;
    for (int i = 0; i < index; i++)
    {
        count += CountSetBits(String[spin][i]);
    }
    long test = (one << bit) - one;
    test &= String[spin][index];
    count += CountSetBits(test);
    return count;
}

//parity
double Determinant::Parity(int orbital, int spin) const
{
    int count = CountSetOrbsTo(orbital, spin);
    if (count % 2 == 0)
        return 1.0;
    else
        return -1.0;
}

//write and read
void Determinant::Write(std::string filename)
{
    std::ofstream ofs(filename, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << *this;
    ofs.close();
}
void Determinant::Read(std::string filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> *this;
    ifs.close();
}

//hartree fock determinant, sets lowest indexed nalpha(nbeta) electrons to occupied
void Determinant::HartreeFock()
{
    for (int i = 0; i < nalpha; i++)
    {
        set(i, 0, true);
    }
    for (int i = 0; i < nbeta; i++)
    {
        set(i, 1, true);
    }
}
