#ifndef HAMILTONIAN_HEADER_H
#define HAMILTONIAN_HEADER_H
#include "Determinant.h"
#include "CIVector.h"
#include "Integrals.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense>



class Hamiltonian
{
    private:
    //boost serialization
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & core_e & norb & nelec & nalpha & nbeta & sz;
        ar & I1 & I2;
        ar & irrep;
    }

    double core_e;
    int norb, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;
    Integral::OneElectron I1;
    Integral::TwoElectron I2;

    public:
    Hamiltonian()
    {
        ReadFCIDUMP("FCIDUMP", I1, I2, core_e, norb, nelec, nalpha, nbeta, sz, irrep);
        InitDetVars(sz, nelec, norb);
    }

    BuildMatrix(const CIVector &V, Eigen::MatrixXd &H)
    {
        H.setZero(V.size, V.size);
        int i = 0, j = 0;
        for (auto it = V.begin(); it != V.end(); ++it)
        {
            for (auto it1 = V.begin(); it1 != V.end(); ++it1)
            {
                H(i, j) = (*this)(it->second, it1->second);
            }
        }
    }

    double operator*(const Determinant &LHS, const Determinant &RHS)
    {
        double H = 0.0;
        //int nelec = Determinant::nelec();
        int n = NumDiffOrbs(LHS, RHS);
        int i = 0, j = 0, a = 0, b = 0;
        if (n == 0)
        {

            std::array<std::vector<int>, 2> open, closed;
            LHS.OpenClosed(open, closed);
            assert(Determinant::nelec == closed.size());
            for (int i = 0; i < closed[0].size(); i++)
            {
                H += I1(closed[0], closed[0]);
                for (int j = 0; j < closed[0]; )
            }
        }
        else if (n == 1)
        {
            OneDiffOrbIndices(LHS, RHS, i, a);
            H += I1(i, a);

        }
        else if (n == 2)
        {
            TwoDiffOrbIndices(LHS, RHS, i, j, a, b);
        }
    }

};
#endif
