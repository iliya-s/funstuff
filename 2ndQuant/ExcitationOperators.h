#ifndef EXCITATION_OPERATORS_HEADER
#define EXCITATION_OPERATORS_HEADER
#include "Determinant.h"
#include "FundamentalOperators.h"

namespace Operator
{
    //Occupation Number Operator
    class Excitation
    {
        private:
        Annihilation a;
        Creation a_dag;

        public:
        //constructors
        Excitation(int orbital_a, int spin_a, int orbital_i, int spin_i) : a(orbital_i, spin_i), a_dag(orbital_a, spin_a) {}
        Excitation(int spin_orbital_a, int spin_orbital_i) : a(spin_orbital_i), a_dag(spin_orbital_a) {}
        Excitation(const Creation _a_dag, Annihilation _a) : a(_a), a_dag(_a_dag) {}

        //setting spin orbital
        Excitation &operator()(int orbital_a, int spin_a, int orbital_i, int spin_i)
        {
            a(orbital_i, spin_i);
            a_dag(orbital_a, spin_a);
            return *this;
        }
        Excitation &operator()(int spin_orbital_a, int spin_orbital_i)
        {
            a(spin_orbital_i);
            a_dag(spin_orbital_a);
            return *this;
        }

        //adjoint
        Excitation adjoint()
        {
            Annihilation c(a_dag);
            Creation c_dag(a);
            return Excitation(c_dag, c);
        }

        //application onto Determinant from left and right
        friend Determinant &operator*(const Excitation &E, Determinant &D);
        friend Determinant &operator*(Determinant &D, const Excitation &E);
    };
    Determinant &operator*(const Excitation &E, Determinant &D)
    {
        D = E.a * D;
        return E.a_dag * D;
    }
    Determinant &operator*(Determinant &D, const Excitation &E)
    {
        D = D * E.a_dag;
        return D * E.a;
    }
}
#endif
