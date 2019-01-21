#ifndef FUNDAMENTAL_OPERATORS_HEADER
#define FUNDAMENTAL_OPERATORS_HEADER
#include "Determinant.h"

namespace Operator
{
    //fundamental operators - Creation and Annihilation Operators
    class Creation;
    class Annihilation;

    class Creation
    {
        private:
        int index; //spin orbital index of Creation operator

        public:
        //constructors
        Creation();
        Creation(int orbital, int spin);
        Creation(int spin_orbital);
        Creation(const Creation &a_dag);
        friend class Annihilation;
        Creation(const Annihilation &a);

        //setting spin orbital
        Creation &operator()(int orbital, int spin);
        Creation &operator()(int spin_orbital);

        //application onto Determinant from left and right
        friend Determinant operator*(const Creation &a_dag, const Determinant &D);
        friend Determinant operator*(const Determinant &D, const Creation &a_dag);
        //excitation operator
        friend class Excitation;
    };

    class Annihilation
    {
        private:
        int index; //spin orbital index of Annihilation operator

        public:
        //constructors
        Annihilation();
        Annihilation(int orbital, int spin);
        Annihilation(int spin_orbital);
        Annihilation(const Annihilation &a);
        friend class Creation;
        Annihilation(const Creation &a_dag);
        
        //setting spin orbital
        Annihilation &operator()(int orbital, int spin);
        Annihilation &operator()(int spin_orbital);

        //application onto Determinant from left and right
        friend Determinant operator*(const Annihilation &a, const Determinant &D);
        friend Determinant operator*(const Determinant &D, const Annihilation &a);
        //excitation operator
        friend class Excitation;
    };
}
#endif
