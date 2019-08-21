#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER
#include "Determinant.h"
#include "FockVector.h"
#include "Integrals.h"
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

/*
class HashPair
{
    public:
    std::size_t operator()(const std::pair<Determinant,Determinant> &ij) const { return std::hash<std::size_t>(ij.first.key() * ij.second.key()); } //I think this will clash in a good way
};
*/

namespace Operator
{
    class Hamiltonian
    {
        private:
        //boost serialization
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & core_e & norb & nelec & nalpha & nbeta & sz;
            ar & irrep;
            ar & I1 & I2;
        }
    
        double core_e;
        int norb, nelec, nalpha, nbeta, sz;
        std::vector<int> irrep;
        Integral::OneElectron I1;
        Integral::TwoElectron I2;
        //std::unordered_map<std::pair<Determinant, Determinant>, double, HashPair> Hij;  //hashtable to store diagonal values and singly connected values

    
        public:
        Hamiltonian(std::string filename = "FCIDUMP")
        {
            ReadFCIDUMP(filename, I1, I2, core_e, norb, nelec, nalpha, nbeta, sz, irrep);
            InitDetVars(sz, nelec, norb);
        }
    
        double energy(const Determinant &D) const
        {
            double E = core_e;
            std::array<std::vector<int>, 2> open, closed;
            D.OpenClosed(open, closed);
            assert(Determinant::nalpha() == closed[0].size());
            assert(Determinant::nbeta() == closed[1].size());
    
            //one electron operator
            for (int sz = 0; sz < 2; sz++) //like spin
            {
                for (int i = 0; i < closed[sz].size(); i++)
                {
                    E += I1(closed[sz][i], closed[sz][i]);
                }
            }
            //two electron operator
            for (int sz = 0; sz < 2; sz++) //like spin
            {
                for (int i = 0; i < closed[sz].size(); i++)
                {
                    for (int j = 0; j < closed[sz].size(); j++)
                    {
                        E += 0.5 * I2(closed[sz][i], closed[sz][i], closed[sz][j], closed[sz][j]); //direct
                        E -= 0.5 * I2(closed[sz][i], closed[sz][j], closed[sz][j], closed[sz][i]); //exchange
                    }
                }
            }
            for (int a = 0; a < closed[0].size(); a++) //differing spin
            {
                for (int b = 0; b < closed[1].size(); b++)
                {
                    E += I2(closed[0][a], closed[0][a], closed[1][b], closed[1][b]); //direct
                }
            }
            return E;
        }
    
        void diagonal(const FockVector &CI, Eigen::VectorXd &V) const
        {
            V.setZero(CI.size());
            int i = 0;
            for (auto it = CI.begin(); it != CI.end(); ++it)
            {
                V(i) = energy(*it);
                i++;
            }
        }
    
        void matrix(const FockVector &CI, Eigen::MatrixXd &H) const
        {
            H.setZero(CI.size(), CI.size());
            int i = 0;
            for (auto it = CI.begin(); it != CI.end(); ++it)
            {
                int j = 0;
                for (auto it1 = CI.begin(); it1 != CI.end(); ++it1)
                {
                    H(i, j) = (*this)(*it, *it1);
                    j++;
                }
                i++;
            }
        }
    
        double operator()(const Determinant &LHS, const Determinant &RHS) const
        {
            double H = 0.0;
            int n = NumDiffOrbs(LHS, RHS);
            if (n == 0) //LHS == RHS
            {
                H = energy(LHS);
            }
            else if (n == 1) //LHS(i -> a) == RHS
            {
                int i = 0, a = 0;
                OneDiffOrbIndices(LHS, RHS, i, a);
    
                int si = i % 2;
                int sa = a % 2;
                int orbi = i / 2;
                int orba = a / 2;
                double pi = LHS.parity(orbi, si);
                double pa = LHS.parity(orba, sa);
    
                //one electron operator
                if (si == sa) //like spin
                {
                    H += pi * pa * I1(orbi, orba);
                }
                //differing spin integral = 0
    
                //two electron operator
                std::vector<int> open, closed;
                LHS.OpenClosed(open, closed);
                for (int r = 0; r < closed.size(); r++)
                {
                    int sr = closed[r] % 2;  
                    int orbr = closed[r] / 2;
    
                    if (si == sa) //direct term
                    {
                        H += pi * pa * I2(orbi, orba, orbr, orbr);
                    }
                    if (si == sr && sa == sr) //exchange term
                    {
                        H -= pi * pa * I2(orbi, orbr, orbr, orba);
                    }
                }
            }
            else if (n == 2) //LHS(ij -> ab) == RHS
            {
                int i = 0, j = 0, a = 0, b = 0;
                TwoDiffOrbIndices(LHS, RHS, i, j, a, b);
    
                int si = i % 2;
                int sj = j % 2;
                int sa = a % 2;
                int sb = b % 2;
                int orbi = i / 2;
                int orbj = j / 2;
                int orba = a / 2;
                int orbb = b / 2;
                double pi = LHS.parity(orbi, si);
                double pj = LHS.parity(orbj, sj);
                double pa = RHS.parity(orba, sa);
                double pb = RHS.parity(orbb, sb);
    
                //one electron operator = 0
                //two electron operator
                if (si == sa && sj == sb) //direct
                {
                    H += pi * pj * pa * pb * I2(orbi, orba, orbj, orbb);
                }
                if (si == sb && sj == sa) //exchange
                {
                    H -= pi * pj * pa * pb * I2(orbi, orbb, orbj, orba);
                }
            }
            return H;
        }
        
        Eigen::VectorXd multiply(const FockVector &CI, const Eigen::VectorXd &z) const
        {
            assert(CI.size() == z.rows());
            Eigen::VectorXd Hz = Eigen::VectorXd::Zero(z.rows());
            int i = 0;
            for (auto row = CI.begin(); row != CI.end(); ++row)
            {
                int j = 0;
                for (auto col = CI.begin(); col != CI.end(); ++col)
                {
                    Hz(i) += (*this)(*row, *col) * z(j);
                    j++;
                }
                i++;
            }
            return Hz;
        }
    };
}
#endif
