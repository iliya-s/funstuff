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
        OccupationNumber(int orbital, int spin) : index(2 * orbital + spin) { assert(index >= 0 && index < 2 * Determinant::norb()); }
        OccupationNumber(int spin_orbital) : index(spin_orbital) { assert(index >= 0 && index < 2 * Determinant::norb()); }

        //setting spin orbital
        OccupationNumber &operator()(int orbital, int spin)
        {
            index = 2 * orbital + spin;
            assert(index >= 0 && index < 2 * Determinant::norb());
            return *this;
        }
        OccupationNumber &operator()(int spin_orbital)
        {
            index = spin_orbital;
            assert(index >= 0 && index < 2 * Determinant::norb());
            return *this;
        }

        //application onto Determinant from left and right
        friend Determinant operator*(const OccupationNumber &n, Determinant D)
        {
            if (!D(n.index)) { D.zero(); } //if index is not occupied -> return 0
            return D;   //if index is occupied -> return Determinant "multiplied by 1.0"
        }
        friend Determinant operator*(Determinant D, const OccupationNumber &n) { return n * D; }
    };


    //Particle Number Operator
    class ParticleNumber
    {
        public:
        //application onto Determinant from left and right
        friend Determinant operator*(const ParticleNumber &N, Determinant D) { return D.CountSetOrbs() * D; }
        friend Determinant operator*(Determinant D, const ParticleNumber &N) { return N * D; }
    };
}
#endif
