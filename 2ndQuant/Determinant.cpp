#include "Determinant.h"
#include <fstream>
#include <iostream>
#include <string>
#include <functional>

//initialize static variables
void InitDetVars(int spin, int nelec, int norb)
{
    Determinant::Norb = norb;
    Determinant::Len = norb / 64 + 1;
    Determinant::Nalpha = nelec / 2 + spin;
    Determinant::Nbeta = nelec - Determinant::Nalpha;
}

//constructors and deconstructor
Determinant::Determinant()
{
    String[0] = new long[Len];
    String[1] = new long[Len];
    vacuum();
}

Determinant::Determinant(const Determinant &D)
{
    Coeff = D.Coeff;
    String[0] = new long[Len];
    String[1] = new long[Len];
    for (int i = 0; i < Len; i++)
    {
        String[0][i] = D.String[0][i];
        String[1][i] = D.String[1][i];
    }
}

Determinant::Determinant(const std::vector<int> &spinOrbs) //builds determinant from vector<int> of occupied spin orbitals
{
    vacuum();
    for (int i = 0; i < spinOrbs.size(); i++)
    {
        set(spinOrbs[i], true);
    }
}

Determinant::Determinant(const std::vector<int> &alpha, const std::vector<int> &beta) //builds determinant from vector<int> of occupied alpha and beta orbitals
{
    String[0] = new long[Len];
    String[1] = new long[Len];
    vacuum();
    for (int i = 0; i < alpha.size(); i++)
    {
        set(alpha[i], 0, true);

    }
    for (int i = 0; i < beta.size(); i++)
    {
        set(beta[i], 1, true);
    }
}

Determinant::~Determinant()
{
    delete[] String[1];
    delete[] String[0];
    delete[] String;
}

//operators

    //unique key to determinant
int Determinant::key() const
{
    std::string str;
    for (int i = 0; i < Len; i++)
    {
        str += std::to_string(String[0][i]) + std::to_string(String[1][i]);
    }
    return std::hash<std::string>{}(str);
}

    //copy/assignment operator
Determinant &Determinant::operator=(Determinant RHS)
{
    Coeff = RHS.Coeff;
    for (int i = 0; i < Len; i++)
    {
        String[0][i] = RHS.String[0][i];
        String[1][i] = RHS.String[1][i];
    }
    return *this;
}

    //multiplication operator overloaded for overlap = <bra|ket>
double Determinant::operator*(const Determinant &RHS) const
{
    for (int i = 0; i < Len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i])
            return 0.0;
    }
    return Coeff * RHS.Coeff;
}

    //multiplication operator overloaded for constants, constant * |ket>, |ket> * constant
Determinant &Determinant::operator*=(double constant)
{
    Coeff *= constant;
    return *this;
}
Determinant operator*(Determinant D, double constant)
{
    return D *= constant;
}
Determinant operator*(double constant, Determinant D)
{
    return D *= constant;
}

    //division operator overloaded for a constant on the right, Determinant / constant
Determinant &Determinant::operator/=(double constant)
{
    Coeff /= constant;
    return *this;
}
Determinant operator/(Determinant D, double constant)
{
    return D /= constant;
}

    //equivalence comparison
bool Determinant::operator==(const Determinant &RHS) const
{
    if (Coeff != RHS.Coeff)
        return false;
    for (int i = 0; i < Len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i])
            return false;
    }
    return true;
}

    //output stream
std::ostream &operator<<(std::ostream &os, const Determinant &D)
{
    os << D.Coeff << " | ";
    for (int i = 0, norb = Determinant::Norb; i < norb; i++)
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

//bit manipulations for spin orbital getter and setter from: https://en.wikipedia.org/wiki/Bit_manipulation under Bit manipulation in the C programming language
//spin orbital getter
bool Determinant::operator()(int orbital, int spin) const
{
    assert(orbital >= 0 && orbital < norb());
    assert(spin == 0 || spin == 1);
    long index = orbital / 64;
    long bit = orbital % 64;
    long test = String[spin][index] & (1 << bit);
    if (test == 0)
        return false;
    else
        return true;
}

bool Determinant::operator()(int spin_orbital) const
{
    int orbital = spin_orbital / 2;
    int spin = spin_orbital % 2;
    return (*this)(orbital, spin);
}

//spin orbital setter
void Determinant::set(int orbital, int spin, bool occupancy)
{
    assert(orbital >= 0 && orbital < norb());
    assert(spin == 0 || spin == 1);
    long index = orbital / 64;
    long bit = orbital % 64;
    if (occupancy)
        String[spin][index] |= (1 << bit);
    else
        String[spin][index] &= ~(1 << bit);
}

void Determinant::set(int spin_orbital, bool occupancy)
{
    int orbital = spin_orbital / 2;
    int spin = spin_orbital % 2;
    set(orbital, spin, occupancy);
}

//coefficient getter and setter
double Determinant::coeff() const
{
    return Coeff;
}
double &Determinant::coeff()
{
    return Coeff;
}
void Determinant::coeff(double coefficient)
{
    Coeff = coefficient;
}

//counts occupied spin orbitals in Determinant
int Determinant::CountSetOrbs() const
{
    int count = 0;
    for (int i = 0; i < Len; i++)
    {
        count += CountSetBits(String[0][i]) + CountSetBits(String[1][i]);
    }
    return count;
}

//Counts occupied spin orbitals up to but not including specified orbital
int Determinant::CountSetOrbsTo(int orbital, int spin) const
{
    assert(orbital >= 0 && orbital < norb());
    assert(spin == 0 || spin == 1);
    long index = orbital / 64;
    long bit = orbital % 64;
    int count = 0;
    for (int i = 0; i < index; i++)
    {
        count += CountSetBits(String[0][i]) + CountSetBits(String[1][i]);
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
double Determinant::parity(int orbital, int spin) const
{
    int count = CountSetOrbsTo(orbital, spin);
    if (count % 2 == 0)
        return 1.0;
    else
        return -1.0;
}

double Determinant::parity(int spin_orbital) const
{
    int count = CountSetOrbsTo(spin_orbital);
    if (count % 2 == 0)
        return 1.0;
    else
        return -1.0;
}
    
//get open and closed orbitals
void Determinant::OpenClosed(std::vector<int> &open, std::vector<int> &closed) const
{
    for (int i = 0; i < 2 * norb(); i++)
    {
        if ((*this)(i)) { closed.push_back(i); }
        else { open.push_back(i); }
    }
}
void Determinant::OpenClosed(std::array<std::vector<int>, 2> &open, std::array<std::vector<int>, 2> &closed) const
{
    for (int i = 0; i < norb(); i++)
    {
        //alpha
        if ((*this)(i, 0)) { closed[0].push_back(i); }
        else { open[0].push_back(i); }

        //beta
        if ((*this)(i, 1)) { closed[1].push_back(i); }
        else { open[1].push_back(i); }
    }
}

//write and read
void Determinant::write(std::string filename)
{
    std::ofstream ofs(filename, std::ios::binary);
    boost::archive::binary_oarchive save(ofs);
    save << *this;
    ofs.close();
}
void Determinant::read(std::string filename)
{
    std::ifstream ifs(filename, std::ios::binary);
    boost::archive::binary_iarchive load(ifs);
    load >> *this;
    ifs.close();
}

//hartree fock determinant, sets lowest indexed nelec electrons to occupied
void Determinant::HartreeFock()
{
    vacuum();
    int nelec = nalpha() + nbeta();
    for (int i = 0; i < nelec; i++)
    {
        set(i, true);
    }
}

//Vacuum vector, sets all orbitals to unoccupied and the coefficient to 1.0
void Determinant::vacuum()
{
    Coeff = 1.0;
    for (int i = 0; i < Len; i++)
    {
        String[0][i] = 0;
        String[1][i] = 0;
    }
}

//zero, sets all orbitals to unoccupied and sets coefficient to 0.0
void Determinant::zero()
{
    Coeff = 0.0;
    for (int i = 0; i < Len; i++)
    {
        String[0][i] = 0;
        String[1][i] = 0;
    }
}

void GenerateCombinations(int n, int k, std::vector<std::vector<int>> &combinations)
{
    combinations.clear();
    //std::cout << n;
    //std::cout << k << std::endl;

    std::vector<bool> v(n);
    std::fill(v.begin(), v.begin() + k, true);
    do
    {
        std::vector<int> comb;
        for (int i = 0; i < n; ++i)
        {
            if (v[i])
            {
                //std::cout << (i + 1) << " ";
                comb.push_back(i);
            }
        }
        combinations.push_back(comb);
        //std::cout << "\n";
    } while (std::prev_permutation(v.begin(), v.end()));
}

//calculates the number of different occupied orbitals between two determinants
int NumDiffOrbs(const Determinant &LHS, const Determinant &RHS)
{
    int count = 0;
    for (int i = 0; i < Determinant::Len; i++)
    {
        count += CountSetBits(LHS.String[0][i] ^ RHS.String[0][i]);
        count += CountSetBits(LHS.String[1][i] ^ RHS.String[1][i]);
    }
    return count / 2;
}

//function to find differing orbital indices inspired from https://lemire.me/blog/2018/02/21/iterating-over-set-bits-quickly/

//finds spin orbital indices of differing occupied orbital for two determinants with 1 differing orbital
//  i corresponds to the spin orbital occupied in LHS, a corresponds to the spin orbital occupied in RHS
void OneDiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &a)
{
    int indices[2];
    for (int l = 0; l < RHS.Len; l++)
    {
        long alpha = LHS.String[0][l] ^ RHS.String[0][l];
        long beta = LHS.String[1][l] ^ RHS.String[1][l];
        int pos = 0;
        while (alpha != 0 && beta != 0)
        {
            indices[pos] = 2 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
            indices[pos] = 2 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
        while (alpha != 0)
        {
            indices[pos] = 2 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
        }
        while (beta != 0)
        {
            indices[pos] = 2 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
    }
    if (LHS(indices[0])) //LHS has ith orbital occupied
    {
        i = indices[0];
        a = indices[1];
    }
    else
    {
        i = indices[1];
        a = indices[0];
    }   
}

//finds spin orbital indices of differing occupied orbitals for two determinants with 2 differing orbitals occupied
//  i, j corresponds to the spin orbitals occupied in LHS, a, b corresponds to the spin orbitals occupied in RHS
void TwoDiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &j, int &a, int &b)
{
    int indices[4];
    for (int l = 0; l < LHS.Len; l++)
    {
        long alpha = LHS.String[0][l] ^ RHS.String[0][l];
        long beta = LHS.String[1][l] ^ RHS.String[1][l];
        int pos = 0;
        while (alpha != 0 && beta != 0)
        {
            indices[pos] = 2 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
            indices[pos] = 2 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
        while (alpha != 0)
        {
            indices[pos] = 2 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
        }
        while (beta != 0)
        {
            indices[pos] = 2 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
    }
    if (LHS(indices[0]) && LHS(indices[1]))
    {
        i = std::min(indices[0], indices[1]);
        j = std::max(indices[0], indices[1]);
        a = std::min(indices[2], indices[3]);
        b = std::max(indices[2], indices[3]);
    }
    else if (LHS(indices[0]) && LHS(indices[2]))
    {
        i = std::min(indices[0], indices[2]);
        j = std::max(indices[0], indices[2]);
        a = std::min(indices[1], indices[3]);
        b = std::max(indices[1], indices[3]);
    }
    else if (LHS(indices[0]) && LHS(indices[3]))
    {
        i = std::min(indices[0], indices[3]);
        j = std::max(indices[0], indices[3]);
        a = std::min(indices[1], indices[2]);
        b = std::max(indices[1], indices[2]);
    }
    else if (LHS(indices[1]) && LHS(indices[2]))
    {
        i = std::min(indices[1], indices[2]);
        j = std::max(indices[1], indices[2]);
        a = std::min(indices[0], indices[3]);
        b = std::max(indices[0], indices[3]);
    }
    else if (LHS(indices[1]) && LHS(indices[3]))
    {
        i = std::min(indices[1], indices[3]);
        j = std::max(indices[1], indices[3]);
        a = std::min(indices[2], indices[0]);
        b = std::max(indices[2], indices[0]);
    }
    else if (LHS(indices[2]) && LHS(indices[3]))
    {
        i = std::min(indices[2], indices[3]);
        j = std::max(indices[2], indices[3]);
        a = std::min(indices[1], indices[0]);
        b = std::max(indices[1], indices[0]);
    }
}
