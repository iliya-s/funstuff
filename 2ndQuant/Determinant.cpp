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
    for (int i = 0; i < len; i++)
    {
        String[0].push_back(0);
        String[1].push_back(0);
    }
}

Determinant::Determinant(const Determinant &D)
{
    coef = D.coef;
    String = D.String;
}

//operators
//  copy/assignment operator
Determinant &Determinant::operator=(const Determinant &RHS)
{
    coef = RHS.coef;
    String = RHS.String;
    return *this;
}

//  multiplication operator overloaded for overlap = <bra|ket>
double Determinant::operator*(const Determinant &RHS) const
{
    for (int i = 0; i < len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i])
            return 0.0;
    }
    return coef * RHS.coef;
}

//  multiplication operator overloaded for constants
Determinant &Determinant::operator*=(double constant)
{
    coef *= constant;
    return *this;
}
Determinant &operator*(Determinant &D, double constant)
{
    return D *= constant;
}
Determinant &operator*(double constant, Determinant &D)
{
    return D *= constant;
}

//  equivalence comparison
bool Determinant::operator==(const Determinant &RHS) const
{
    if (coef != RHS.coef)
        return false;
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
    os << D.coef << " | ";
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
    os << ">";
    return os;
}

//bit manipulations for getter and setter from: https://en.wikipedia.org/wiki/Bit_manipulation under Bit manipulation in the C programming language
//getter
inline bool Determinant::operator()(int orbital, int spin) const
{
    long index = orbital / 64;
    long bit = orbital % 64;
    long test = String[spin][index] & (1 << bit);
    if (test == 0)
        return false;
    else
        return true;
}

inline bool Determinant::operator()(int spin_orbital) const
{
    int orbital = spin_orbital / 2;
    int spin = spin_orbital % 2;
    return (*this)(orbital, spin);
}

//setter
inline void Determinant::set(int orbital, int spin, bool occupancy)
{
    long index = orbital / 64;
    long bit = orbital % 64;
    if (occupancy)
        String[spin][index] |= (1 << bit);
    else
        String[spin][index] &= ~(1 << bit);
}

inline void Determinant::set(int spin_orbital, bool occupancy)
{
    int orbital = spin_orbital / 2;
    int spin = spin_orbital % 2;
    set(orbital, spin, occupancy);
}

//Counts occupied spin orbitals up to but not including specified orbital
int Determinant::CountSetOrbsTo(int orbital, int spin) const
{
    long index = orbital / 64;
    long bit = orbital % 64;
    int count = 0;
    for (int i = 0; i < index; i++)
    {
        count += CountSetBits(String[0][i]);
        count += CountSetBits(String[1][i]);
    }
    long alpha, beta;
    if (spin == 0)  //alpha orbital
    {
        alpha = (1 << bit) - 1;
        beta = alpha;
    }
    else    //beta orbital
    {
        alpha = (1 << (bit + 1)) - 1;
        beta = (1 << bit) - 1;
    }
    alpha &= String[0][index];
    beta &= String[1][index];
    count += CountSetBits(alpha) + CountSetBits(beta);
    return count;
}

int Determinant::CountSetOrbsTo(int spin_orbital) const
{
    int orbital = spin_orbital / 2;
    int spin = spin_orbital % 2;
    return CountSetOrbsTo(orbital, spin);
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

double Determinant::Parity(int spin_orbital) const
{
    int count = CountSetOrbsTo(spin_orbital);
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

//hartree fock determinant, sets lowest indexed nelec electrons to occupied
void Determinant::HartreeFock()
{
    int nelec = nalpha + nbeta;
    for (int i = 0; i < nelec; i++)
    {
        set(i, true);
    }
}
