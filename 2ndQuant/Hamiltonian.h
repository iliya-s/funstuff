#ifndef HAMILTONIAN_HEADER_H
#define HAMILTONIAN_HEADER_H
#include "Determinant.h"

//calculates the number of different occupied orbitals between two determinants
int NumDiffOrbs(const Determinant &LHS, const Determinant &RHS);

//finds spin orbital indices of differing occupied orbital for two determinants with 1 differing orbital
//  i corresponds to the spin orbital occupied in LHS, a corresponds to the spin orbital occupied in RHS
void OneDiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &a);

//finds spin orbital indices of differing occupied orbitals for two determinants with 2 differing orbitals occupied
//  i, j corresponds to the spin orbitals occupied in LHS, a, b corresponds to the spin orbitals occupied in RHS
void TwoDiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &j, int &a, int &b);

/*
class Hamiltonian
{
    private:
    //boost serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
    }

    public:
    Hamiltonian(

};
*/
/*
void DiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &j, int &a, int &b)
{
    int n = NumDiffOrbs(LHS, RHS);
    if (n == 1)
    {
        OneDiffOrbIndices(LHS, RHS, i, a);
        j = 0;
        b = 0;
    }
    else if (n == 2)
    {
        TwoDiffOrbIndices(LHS, RHS, i, j, a, b);
    }
    else
    {
        i = 0;
        j = 0;
        a = 0;
        b = 0;
    }
}
*/
#endif
