#ifndef DETERMINANT_HEADER
#define DETERMINANT_HEADER
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
void InitDetVars(int spin, int nelec, int norb);

class Determinant
{
    private:
    //boost serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & Coeff;
        for (int i = 0; i < Len; i++) { ar & String[0][i] & String[1][i]; }
    }

    //initialize determinant static variables
    friend void InitDetVars(int spin, int nelec, int norb);

    //member variables
    static int Norb, Nalpha, Nbeta, Len;
    mutable double Coeff = 1.0;
    long **String = new long *[2];

    public:
    //return static values
    static int norb() { return Norb; }
    static int nalpha() { return Nalpha; }
    static int nbeta() { return Nbeta; }
    //constructors and deconstructor
    Determinant();
    Determinant(const Determinant &D);
    Determinant(const std::vector<int> &spinOrbs); //builds determinant from vector<int> of occupied spin orbitals
    Determinant(const std::vector<int> &alpha, const std::vector<int> &beta); //builds determinant from vector<int> of occupied alpha and beta orbitals
    ~Determinant();
    
    //operators
    Determinant &operator=(Determinant RHS); //copy/assignment operator
    //multiplication operator overloaded for <bra|ket>, constant * |ket>, |ket> * constant
    double operator*(const Determinant &RHS) const;
    Determinant &operator*=(double constant);
    friend Determinant operator*(Determinant D, double constant);
    friend Determinant operator*(double constant, Determinant D);
    //division operator overloaded for Determinant / constant
    Determinant &operator/=(double constant);
    friend Determinant operator/(Determinant D, double constant);
    bool operator==(const Determinant &RHS) const; //equivalence comparison
    bool operator<(const Determinant &RHS) const; //less than comparison
    friend std::ostream &operator<<(std::ostream &os, const Determinant &D); //output stream

    //functions
    std::size_t key() const; //unique key to determinant
    //spin orbital getter and setter
    bool operator()(int orbital, int spin) const;
    bool operator()(int spin_orbital) const;
    void set(int orbital, int spin, bool occupancy);
    void set(int spin_orbital, bool occupancy);
    //coefficient getter and setter;
    double coeff() const;
    void coeff(double coefficient);

    //counts occupied spin orbitals in Determinant
    int CountSetOrbs() const;
    //counts occupied spin orbitals up to but not including specified orbital
    int CountSetOrbsTo(int orbital, int spin) const;
    int CountSetOrbsTo(int spin_orbital) const;
    //parity
    double parity(int orbital, int spin) const;
    double parity(int spin_orbital) const;
    //get open and closed orbitals
    void OpenClosed(std::vector<int> &open, std::vector<int> &closed) const;
    void OpenClosed(std::array<std::vector<int>, 2> &open, std::array<std::vector<int>, 2> &closed) const;

    //write and read
    void write(std::string filename = "Determinant.bkp");
    void read(std::string filename = "Determinant.bkp");
    
    void HartreeFock(); //hartree fock determinant, sets lowest indexed nalpha and nbeta orbitals to occupied
    void vacuum(); //vacuum vector, sets all orbitals to unoccupied and the coefficient to 1.0
    void zero(); //zero, sets all orbitals to unoccupied and sets coefficient to 0.0
    //generates all connected determinants from a given det, ie d + single excitations + double excitations
    int numConnected() const;
    void connected(std::vector<Determinant> &dets) const;
    //generates all singly connected determinants from a given det, ie single excitations
    int numSinglyConnected() const;
    void singlyConnected(std::vector<Determinant> &dets) const;
    //generates all doubly connected determinants from a given det, ie double excitations
    int numDoublyConnected() const;
    void doublyConnected(std::vector<Determinant> &dets) const;

    //friend functions for hamiltonian overlap
    friend class Hamiltonian;
    friend class FockVector; //friend object Fock vector
    friend int NumDiffOrbs(const Determinant &LHS, const Determinant &RHS); //calculates the number of different occupied orbitals between two determinants
    friend void OneDiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &a); //finds spin orbital indices of differing occupied orbital for two determinants with 1 differing orbital
    friend void TwoDiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &j, int &a, int &b); //finds spin orbital indices of differing occupied orbitals for two determinants with 2 differing orbitals occupied
};

//generates all n choose k combinations of integers and stores them in combinations
void GenerateCombinations(int n, int k, std::vector<std::vector<int>> &combinations);
#endif
