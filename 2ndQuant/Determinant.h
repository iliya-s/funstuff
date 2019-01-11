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
            ar & coef;
            for (int i = 0; i < len; i++)
            {
                ar & String[0][i] & String[1][i];
            }
        }

    //initialize determinant static variables
    friend void InitDetVars(int spin, int nelec, int norbs);

    //member variables
    static int norbs, nalpha, nbeta, len;
    double coef = 1.0;
    std::array<std::vector<long>, 2> String;

    public:
    //constructors and deconstructor
    Determinant();
    Determinant(const Determinant &D);
    
    //operators
    //  copy/assignment operator
    Determinant &operator=(const Determinant &RHS);
    //  multiplication operator overloaded for <bra|ket>, and const * |ket>
    double operator*(const Determinant &RHS) const;
    Determinant &operator*=(double constant);
    //  equivalence comparison
    bool operator==(const Determinant &RHS) const;
    //  output stream
    friend std::ostream &operator<<(std::ostream &os, const Determinant &D);

    //getter and setter
    bool operator()(int orbital, int spin) const;
    bool operator()(int spin_orbital) const;
    void set(int orbital, int spin, bool occupancy);
    void set(int spin_orbital, bool occupancy);

    //counts occupied spin orbitals up to but not including specified orbital
    int CountSetOrbsTo(int orbital, int spin) const;
    int CountSetOrbsTo(int spin_orbital) const;
    //parity
    double Parity(int orbital, int spin) const;
    double Parity(int spin_orbital) const;

    //write and read
    void Write(std::string filename = "Determinant.bkp");
    void Read(std::string filename = "Determinant.bkp");
    
    //hartree fock determinant, sets lowest indexed nelec spin orbitals to occupied
    void HartreeFock();
};
Determinant &operator*(Determinant &D, double constant);
Determinant &operator*(double constant, Determinant &D);
#endif
