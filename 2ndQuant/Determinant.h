#ifndef DETERMINANT_HEADER_H
#define DETERMINANT_HEADER_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#define MAXLONG 0xFFFFFFFFFFFFFFFF //all 64 bits set

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
    Determinant &operator=(const Determinant &RHS);
        //multiplication operator overloaded for <bra|ket>, and constant * |ket>
    double operator*(const Determinant &RHS) const;
    Determinant &operator*=(double constant);
    friend Determinant &operator*(Determinant &D, double constant);
    friend Determinant &operator*(double constant, Determinant &D);
        //equivalence comparison
    bool operator==(const Determinant &RHS) const;
        //output stream
    friend std::ostream &operator<<(std::ostream &os, const Determinant &D);

    //getter and setter
    bool operator()(int orbital, int spin) const;
    bool operator()(int spin_orbital) const;
    void set(int orbital, int spin, bool occupancy);
    void set(int spin_orbital, bool occupancy);
        //getter for coefficient
    double coeff();

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
    //zero, sets all orbitals to occupied and sets coefficient to 0.0
    void zero();
};
#endif
