#ifndef NUMBER_OPERATORS_HEADER
#define NUMBER_OPERATORS_HEADER
#include "Determinant.h"
#include "FundamentalOperators.h"

namespace Operator
{
    //Occupation Number Operator
    class OccupationNumber
    {
        private:
        int index; //spin orbital index of Number operator

        public:
        //constructors
        OccupationNumber() : index(-1) {}
        OccupationNumber(int orbital, int spin) : index(2 * orbital + spin) {}
        OccupationNumber(int spin_orbital) : index(spin_orbital) {}

        //setting spin orbital
        OccupationNumber &operator()(int orbital, int spin)
        {
            index = 2 * orbital + spin;
            return *this;
        }
        OccupationNumber &operator()(int spin_orbital)
        {
            index = spin_orbital;
            return *this;
        }

        //application onto Determinant from left and right
        friend Determinant operator*(const OccupationNumber &n, const Determinant &D)
        {
            Determinant Dcopy(D);
            if (!Dcopy(n.index)) //if index is not occupied -> return 0
                Dcopy.zero();
            return Dcopy;   //if index is occupied -> return Determinant "multiplied by 1.0"
        }
        friend Determinant operator*(const Determinant &D, const OccupationNumber &n)
        {
            return n * D;
        }
    };


    //Particle Number Operator
    class ParticleNumber
    {
        public:
        //application onto Determinant from left and right
        friend Determinant operator*(const ParticleNumber &N, const Determinant &D)
        {
            int count = D.CountSetOrbs();
            return count * D;
        }
        friend Determinant operator*(const Determinant &D, const ParticleNumber &N)
        {
            return N * D;
        }
    };
}
#endif
