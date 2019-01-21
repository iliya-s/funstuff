#include "Determinant.h"
#include "FundamentalOperators.h"

namespace Operator
{
    //Creation Operator
        //constructors
    Creation::Creation() : index(-1) {}
    Creation::Creation(int orbital, int spin) : index(2 * orbital + spin) {}
    Creation::Creation(int spin_orbital) : index(spin_orbital) {}
    Creation::Creation(const Creation &a_dag) : index(a_dag.index) {}
    Creation::Creation(const Annihilation &a) : index(a.index) {}

        //setting spin orbital
    Creation &Creation::operator()(int orbital, int spin)
    {
        index = 2 * orbital + spin;
        return *this;
    }
    Creation &Creation::operator()(int spin_orbital)
    {
        index = spin_orbital;
        return *this;
    }
    
        //application onto Determinant from left and right
    Determinant operator*(const Creation &a_dag, const Determinant &D)
    {
        Determinant Dcopy(D);
        if (Dcopy(a_dag.index))    //if spin orbital is occupied in D -> return 0
            Dcopy.zero();
        else    //if spin orbital is unocuupied in D -> return D with spin orbital occupied and parity multiplied to coefficient
        {
            Dcopy *= Dcopy.parity(a_dag.index);
            Dcopy.set(a_dag.index, true);
        }
        return Dcopy;
    } 
    Determinant operator*(const Determinant &D, const Creation &a_dag)   //multiplcation on the right by a Creation Operator is quivalent to multiplication on the left by an Annihilation Operator
    {
        Annihilation a(a_dag);
        return a * D;
    }


    //Annihilation Operator
        //constructors
    Annihilation::Annihilation() : index(-1) {}
    Annihilation::Annihilation(int orbital, int spin) : index(2 * orbital + spin) {}
    Annihilation::Annihilation(int spin_orbital) : index(spin_orbital) {}
    Annihilation::Annihilation(const Annihilation &a) : index(a.index) {}
    Annihilation::Annihilation(const Creation &a_dag) : index(a_dag.index) {}
    
        //setting spin orbital
    Annihilation &Annihilation::operator()(int orbital, int spin)
    {
        index = 2 * orbital + spin;
        return *this;
    }
    Annihilation &Annihilation::operator()(int spin_orbital)
    {
        index = spin_orbital;
        return *this;
    }
    
        //application onto Determinant from left and right
    Determinant operator*(const Annihilation &a, const Determinant &D)
    {
        Determinant Dcopy(D);
        if (Dcopy(a.index)) //if spin orbital is occupied in D -> return D with spin orbital unoccupied and parity mutliplied to coefficient
        {
            Dcopy *= Dcopy.parity(a.index);
            Dcopy.set(a.index, false);
        }
        else    //if spin orbital is unoccupied in D -> return 0
            Dcopy.zero();
        return Dcopy;
    }

    Determinant operator*(const Determinant &D, const Annihilation &a)   //multiplcation on the right by an Annihilation Operator is quivalent to multiplication on the left by a Creation Operator
    {
        Creation a_dag(a);
        return a_dag * D;
    }
}
