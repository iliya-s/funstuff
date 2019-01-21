#ifndef DETERMINANT_HEADER_H
#define DETERMINANT_HEADER_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#define MAXLONG 0xFFFFFFFFFFFFFFFF //all 64 bits set

//counts number of set bits in integer type
//  uses SWAR algorithim from: https://www.playingwithpointers.com/blog/swar.html
inline int CountSetBits(long n)
{
    n = n - ((n >> 1) & 0x5555555555555555);
    n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
    n = (n + (n >> 4) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101;
    return n >> 56;
}

//initialize determinant static variables
void InitDetVars(int spin, int nelec, int norbs);

class Determinant
{
    private:
    //boost serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & Coeff;
        for (int i = 0; i < len; i++)
        {
            ar & String[0][i] & String[1][i];
        }
    }

    //initialize determinant static variables
    friend void InitDetVars(int spin, int nelec, int norbs);

    //member variables
    static int norbs, nalpha, nbeta, len;
    double Coeff = 1.0;
    long **String = new long *[2];

    public:
    //constructors and deconstructor
    Determinant();
    Determinant(const Determinant &D);
    ~Determinant();
    
    //operators
        //copy/assignment operator
    Determinant &operator=(Determinant RHS);
        //multiplication operator overloaded for <bra|ket>, constant * |ket>, |ket> * constant
    double operator*(const Determinant &RHS) const;
    Determinant &operator*=(double constant);
    friend Determinant operator*(const Determinant &D, double constant);
    friend Determinant operator*(double constant, const Determinant &D);
        //division operator overloaded for Determinant / constant
    Determinant &operator/=(double constant);
    friend Determinant operator/(const Determinant &D, double constant);
        //equivalence comparison
    bool operator==(const Determinant &RHS) const;
        //output stream
    friend std::ostream &operator<<(std::ostream &os, const Determinant &D);

    //spin orbital getter and setter
    bool operator()(int orbital, int spin) const;
    bool operator()(int spin_orbital) const;
    void set(int orbital, int spin, bool occupancy);
    void set(int spin_orbital, bool occupancy);
    //coefficient getter and setter;
    double coeff();
    void coeff(double coefficient);

    //counts occupied spin orbitals in Determinant
    int CountSetOrbs() const;
    //counts occupied spin orbitals up to but not including specified orbital
    int CountSetOrbsTo(int orbital, int spin) const;
    int CountSetOrbsTo(int spin_orbital) const;
    //parity
    double parity(int orbital, int spin) const;
    double parity(int spin_orbital) const;

    //write and read
    void write(std::string filename = "Determinant.bkp");
    void read(std::string filename = "Determinant.bkp");
    
    //hartree fock determinant, sets lowest indexed nelec spin orbitals to occupied
    void HartreeFock();
    //vacuum vector, sets all orbitals to unoccupied and the coefficient to 1.0
    void vacuum();
    //zero, sets all orbitals to unoccupied and sets coefficient to 0.0
    void zero();

    //friend functions for hamiltonian overlap
    friend int NumDiffOrbs(const Determinant &LHS, const Determinant &RHS);
    friend void DiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &a);
    friend void DiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &j, int &a, int &b);
};
#endif
