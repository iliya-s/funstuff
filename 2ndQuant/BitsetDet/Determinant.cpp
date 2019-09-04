#include "Determinant.h"
#include "FockVector.h"
#include <fstream>
#include <string>
#include <functional>

//initialize static variables
void InitDetVars(int spin, int nelec, int norb)
{
    Determinant::Norb = norb;
    Determinant::Len = norb / 64 + 1;
    Determinant::Nalpha = (nelec + spin) / 2;
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
    String[0] = new long[Len];
    String[1] = new long[Len];
    vacuum();
    for (int i = 0; i < spinOrbs.size(); i++) { set(spinOrbs[i], true); }
}

Determinant::Determinant(const std::vector<int> &alpha, const std::vector<int> &beta) //builds determinant from vector<int> of occupied alpha and beta orbitals
{
    String[0] = new long[Len];
    String[1] = new long[Len];
    vacuum();
    for (int i = 0; i < alpha.size(); i++) { set(alpha[i], 0, true); }
    for (int i = 0; i < beta.size(); i++) { set(beta[i], 1, true); }
}

Determinant::~Determinant()
{
    delete[] String[1];
    delete[] String[0];
    delete[] String;
}

//operators

std::size_t Determinant::key() const //unique key to determinant
{
    std::string alpha, beta;
    for (int i = 0; i < Len; i++)
    {
        alpha += std::to_string(String[0][i]) + std::to_string(i);
        beta += std::to_string(String[1][i]) + std::to_string(i);
    }
    return std::hash<std::string>{}(alpha + beta);
}

Determinant &Determinant::operator=(Determinant RHS) //copy/assignment operator
{
    Coeff = RHS.Coeff;
    for (int i = 0; i < Len; i++)
    {
        String[0][i] = RHS.String[0][i];
        String[1][i] = RHS.String[1][i];
    }
    return *this;
}

double Determinant::operator*(const Determinant &RHS) const //multiplication operator overloaded for overlap = <bra|ket>
{
    for (int i = 0; i < Len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i])
            return 0.0;
    }
    return Coeff * RHS.Coeff;
}

Determinant &Determinant::operator*=(double constant) //multiplication operator overloaded for constants, constant * |ket>, |ket> * constant
{
    Coeff *= constant;
    return *this;
}
Determinant operator*(Determinant D, double constant) { return D *= constant; }
Determinant operator*(double constant, Determinant D) { return D *= constant; }

Determinant &Determinant::operator/=(double constant) //division operator overloaded for a constant on the right, Determinant / constant
{
    Coeff /= constant;
    return *this;
}
Determinant operator/(Determinant D, double constant) { return D /= constant; }

bool Determinant::operator==(const Determinant &RHS) const //equivalence comparison
{
    bool test = true;
    if (Coeff != RHS.Coeff) { test == false; }
    for (int i = 0; i < Len; i++)
    {
        if (String[0][i] != RHS.String[0][i] || String[1][i] != RHS.String[1][i]) { test = false; }
    }
    return test;
}

bool Determinant::operator<(const Determinant &RHS) const //less than comparison
{
    bool test = false;
    for (int i = Len - 1; i >= 0; i--)
    {
        if (String[0][i] < RHS.String[0][i])
        {
            test = true;
            break;
        }
        else if (String[0][i] > RHS.String[0][i]) { break; }
        
        if (String[1][i] < RHS.String[1][i])
        {
            test = true;
            break;
        }
        else if (String[1][i] > RHS.String[1][i]) { break; }
    }
    return test;
}

std::ostream &operator<<(std::ostream &os, const Determinant &D) //output stream
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
    if (test == 0) { return false; }
    else { return true; }
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
    if (occupancy) { String[spin][index] |= (1 << bit);}
    else { String[spin][index] &= ~(1 << bit); }
}

void Determinant::set(int spin_orbital, bool occupancy)
{
    int orbital = spin_orbital / 2;
    int spin = spin_orbital % 2;
    set(orbital, spin, occupancy);
}

//coefficient getter and setter
double Determinant::coeff() const { return Coeff; }
void Determinant::coeff(double coefficient) { Coeff = coefficient; }

//counts occupied spin orbitals in Determinant
int Determinant::CountSetOrbs() const
{
    int count = 0;
    for (int i = 0; i < Len; i++) { count += CountSetBits(String[0][i]) + CountSetBits(String[1][i]); }
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
    for (int i = 0; i < index; i++) { count += CountSetBits(String[0][i]) + CountSetBits(String[1][i]); }
    long alpha, beta;
    if (spin == 0) //alpha orbital
    {
        alpha = (1 << bit) - 1;
        beta = alpha;
    }
    else //beta orbital
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
    if (count % 2 == 0) { return 1.0; }
    else { return -1.0; }
}

double Determinant::parity(int spin_orbital) const
{
    int count = CountSetOrbsTo(spin_orbital);
    if (count % 2 == 0) { return 1.0; }
    else { return -1.0; }
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
void Determinant::write(std::string filename) const
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

//hartree fock determinant, sets lowest indexed nalpha and nbeta electrons to occupied
void Determinant::HartreeFock()
{
    vacuum();
    for (int i = 0; i < nalpha(); i++) { set(i, 0, true); }
    for (int i = 0; i < nbeta(); i++) { set(i, 1, true); }
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

//generates excitations
int Determinant::nSingleExcitations() const
{
    int num = 0;
    std::array<std::vector<int>, 2> open, closed;
    OpenClosed(open, closed);
    for (int sz = 0; sz < 2; sz++) { num += closed[sz].size() * open[sz].size(); }
    return num;
}

void Determinant::singleExcitations(std::vector<Determinant> &dets) const
{
    std::array<std::vector<int>, 2> open, closed;
    OpenClosed(open, closed);

    for (int sz = 0; sz < 2; sz++)
    {
        for (int i = 0; i < closed[sz].size(); i++)
        {
            for (int a = 0; a < open[sz].size(); a++)
            {
                int p = closed[sz][i];
                int q = open[sz][a];
                Determinant dexcite(*this);
                dexcite.set(p, sz, false);
                dexcite.set(q, sz, true); 
                dets.push_back(dexcite);
            }
        }
    }    
}

void Determinant::singleExcitations(std::vector<std::pair<int, int>> &excite) const
{
    std::array<std::vector<int>, 2> open, closed;
    OpenClosed(open, closed);

    for (int sz = 0; sz < 2; sz++)
    {
        for (int i = 0; i < closed[sz].size(); i++)
        {
            for (int a = 0; a < open[sz].size(); a++)
            {
                int p = 2 * closed[sz][i] + sz;
                int q = 2 * open[sz][a] + sz;
                auto pq = std::make_pair(p, q);
                excite.push_back(pq);
            }
        }
    }    
}

int Determinant::nDoubleExcitations() const
{
    int num = 0;
    std::array<std::vector<int>, 2> open, closed;
    OpenClosed(open, closed);

    for (int sz = 0; sz < 2; sz++) //like spin
    {
        int pickclosed = closed[sz].size() * (closed[sz].size() - 1) / 2;
        int pickopen = open[sz].size() * (open[sz].size() - 1) / 2;
        num +=  pickclosed * pickopen;
    }
    //differing spin
    int pickclosed = closed[0].size() * closed[1].size();
    int pickopen = open[0].size() * open[1].size();    
    num += pickclosed * pickopen;
    return num;
}

void Determinant::doubleExcitations(std::vector<Determinant> &dets) const
{
    std::array<std::vector<int>, 2> open, closed;
    OpenClosed(open, closed);

    for (int sz = 0; sz < 2; sz++) //like spin
    {
        for (int i = 0; i < closed[sz].size(); i++)
        {
            for (int j = i + 1; j < closed[sz].size(); j++)
            {
                for (int a = 0; a < open[sz].size(); a++)
                {
                    for (int b = a + 1; b < open[sz].size(); b++)
                    {
                        int p = closed[sz][i];
                        int q = closed[sz][j];
                        int r = open[sz][a];
                        int s = open[sz][b];
                        Determinant dexcite(*this);
                        dexcite.set(p, sz, false);
                        dexcite.set(q, sz, false); 
                        dexcite.set(r, sz, true);
                        dexcite.set(s, sz, true); 
                        dets.push_back(dexcite);
                    }
                }
            }
        }
    }
    for (int i = 0; i < closed[0].size(); i++) //differing spin
    {
        for (int j = 0; j < closed[1].size(); j++)
        {
            for (int a = 0; a < open[0].size(); a++)
            {
                for (int b = 0; b < open[1].size(); b++)
                {
                    int p = closed[0][i];
                    int q = closed[1][j];
                    int r = open[0][a];
                    int s = open[1][b];
                    Determinant dexcite(*this);
                    dexcite.set(p, 0, false);
                    dexcite.set(q, 1, false); 
                    dexcite.set(r, 0, true);
                    dexcite.set(s, 1, true); 
                    dets.push_back(dexcite);
                }
            }
        }
    }
}
    
void Determinant::doubleExcitations(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> &excite) const
{
    std::array<std::vector<int>, 2> open, closed;
    OpenClosed(open, closed);

    for (int sz = 0; sz < 2; sz++) //like spin
    {
        for (int i = 0; i < closed[sz].size(); i++)
        {
            for (int j = i + 1; j < closed[sz].size(); j++)
            {
                for (int a = 0; a < open[sz].size(); a++)
                {
                    for (int b = a + 1; b < open[sz].size(); b++)
                    {
                        int p = 2 * closed[sz][i] + sz;
                        int q = 2 * closed[sz][j] + sz;
                        int r = 2 * open[sz][a] + sz;
                        int s = 2 * open[sz][b] + sz;
                        auto pq = std::make_pair(p, q);
                        auto rs = std::make_pair(r, s);
                        excite.push_back(std::make_pair(pq, rs));
                    }
                }
            }
        }
    }
    for (int i = 0; i < closed[0].size(); i++) //differing spin
    {
        for (int j = 0; j < closed[1].size(); j++)
        {
            for (int a = 0; a < open[0].size(); a++)
            {
                for (int b = 0; b < open[1].size(); b++)
                {
                    int p = 2 * closed[0][i] + 0;
                    int q = 2 * closed[1][j] + 1;
                    int r = 2 * open[0][a] + 0;
                    int s = 2 * open[1][b] + 1;
                    auto pq = std::make_pair(p, q);
                    auto rs = std::make_pair(r, s);
                    excite.push_back(std::make_pair(pq, rs));
                }
            }
        }
    }
}

int Determinant::nExcitations() const { return nSingleExcitations() + nDoubleExcitations(); }

void Determinant::excitations(std::vector<Determinant> &dets) const
{
    singleExcitations(dets);
    doubleExcitations(dets);
}

void Determinant::screenedSingleExcitations(const Integral::HeatBath::OneElectron &HBI1, std::vector<Determinant> &dets, double epsilon) const
{
    std::vector<int> open, closed;
    OpenClosed(open, closed);
    for (int i = 0; i < closed.size(); i++)
    {
        for (auto it = HBI1.begin(closed[i]); it != HBI1.end(closed[i]); ++it)
        {
            double val = it->first;
            int a = it->second;
            if ((*this)(a)) { continue; }
            if (val < epsilon) { break; }
            else
            {
                Determinant dexcite(*this);
                dexcite.set(closed[i], false);
                dexcite.set(a, true);
                dets.push_back(dexcite);
            }
        }
    }
}

void Determinant::screenedDoubleExcitations(const Integral::HeatBath::TwoElectron &HBI2, std::vector<Determinant> &dets, double epsilon) const
{
    std::vector<int> open, closed;
    OpenClosed(open, closed);
    for (int i = 0; i < closed.size(); i++)
    {
        for (int j = i + 1; j < closed.size(); j++)
        {
            auto ij = std::make_pair(closed[i], closed[j]);
            for (auto it = HBI2.begin(ij); it != HBI2.end(ij); ++it)
            {
                double val = it->first;
                auto ab = it->second;
                int a = ab.first;
                int b = ab.second;
                if ((*this)(a) || (*this)(b)) { continue; }
                if (val < epsilon) { break; }
                else
                {
                    Determinant dexcite(*this);
                    dexcite.set(closed[i], false);
                    dexcite.set(closed[j], false);
                    dexcite.set(a, true);
                    dexcite.set(b, true);
                    dets.push_back(dexcite);
                }
            }
        }
    }

}

void Determinant::screenedExcitations(const Integral::HeatBath::OneElectron &HBI1, const Integral::HeatBath::TwoElectron &HBI2, std::vector<Determinant> &dets, double epsilon) const
{
    screenedSingleExcitations(HBI1, dets, epsilon);
    screenedDoubleExcitations(HBI2, dets, epsilon);
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
std::pair<int, int> OneDiffOrbIndices(const Determinant &LHS, const Determinant &RHS)
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
    int i, a;
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
    return std::make_pair(i, a);
}

//finds spin orbital indices of differing occupied orbitals for two determinants with 2 differing orbitals occupied
//  i, j corresponds to the spin orbitals occupied in LHS, a, b corresponds to the spin orbitals occupied in RHS
std::pair<std::pair<int, int>, std::pair<int, int>> TwoDiffOrbIndices(const Determinant &LHS, const Determinant &RHS)
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
    int i, j, a, b;
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
    return std::pair<std::pair<int, int>, std::pair<int, int>>(std::make_pair(i, j), std::make_pair(a, b));
}

void GenerateCombinations(int n, int k, std::vector<std::vector<int>> &combinations)
{
    combinations.clear();

    std::vector<bool> v(n);
    std::fill(v.begin(), v.begin() + k, true);
    do
    {
        std::vector<int> comb;
        for (int i = 0; i < n; ++i) { if (v[i]) { comb.push_back(i); } }
        combinations.push_back(comb);
    } while (std::prev_permutation(v.begin(), v.end()));
}
