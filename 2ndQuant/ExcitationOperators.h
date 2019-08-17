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
        Excitation() {}
        Excitation(int orbital_a, int spin_a, int orbital_i, int spin_i) : a(2 * orbital_i + spin_i), a_dag(2 * orbital_a + spin_a) {}
        Excitation(int spin_orbital_a, int spin_orbital_i) : a(spin_orbital_i), a_dag(spin_orbital_a) {}
        Excitation(const Creation &_a_dag, const Annihilation &_a) : a(_a), a_dag(_a_dag) {}

        //setting spin orbital
        Excitation &operator()(int orbital_a, int spin_a, int orbital_i, int spin_i)
        {
            a(2 * orbital_i + spin_i);
            a_dag(2 * orbital_a + spin_a);
            return *this;
        }
        Excitation &operator()(int spin_orbital_a, int spin_orbital_i)
        {
            a(spin_orbital_i);
            a_dag(spin_orbital_a);
            return *this;
        }

        //adjoint
        Excitation adjoint() const
        {
            Annihilation(a_dag);
            Creation(a);
            return *this;
        }

        //application onto Determinant from left and right
        friend Determinant operator*(const Excitation &E, const Determinant &D)
        {
            Determinant Dcopy(E.a * D);
            return E.a_dag * Dcopy;
        }
        friend Determinant operator*(const Determinant &D, const Excitation &E)
        {
            return E.adjoint() * D;
        }
    };
}
#endif
