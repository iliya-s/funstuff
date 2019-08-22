#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER
#include "Integrals.h"
#include "Determinant.h"
#include "FockVector.h"
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

struct HashKey { std::size_t operator()(const std::size_t &key) const { return key; } };

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
            ar & Hii & Hia;
        }
    
        double core_e;
        int norb, nelec, nalpha, nbeta, sz;
        std::vector<int> irrep;
        Integral::OneElectron I1;
        Integral::TwoElectron I2;
        mutable std::unordered_map<std::size_t, double, HashKey> Hii;  //hashtable to store diagonal values
        mutable std::unordered_map<std::size_t, double, HashKey> Hia;  //hashtable to store singly connected values
    
        public:
        Hamiltonian(std::string filename = "FCIDUMP")
        {
            ReadFCIDUMP(filename, I1, I2, core_e, norb, nelec, nalpha, nbeta, sz, irrep);
            InitDetVars(sz, nelec, norb);
        }

        void clear() const { Hii.clear(); Hia.clear(); } //clears hashtables 
        std::size_t size() const { return Hii.size() + Hia.size(); } //returns size of hashtables
    
        double energy(const Determinant &D) const
        {
            double E = 0.0;
            //if element has been calculated, extract from hashtable
            auto key = D.key();
            if (Hii.count(key) == 1) { E = Hii.at(key); }
            else
            {
                E += core_e;
                std::array<std::vector<int>, 2> open, closed;
                D.OpenClosed(open, closed);
    
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

                //add newly calculated element to hashtable
                Hii[key] = E; 
            }
            return E;
        }

        double operator()(const Determinant &LHS, const Determinant &RHS) const
        {
            double H = 0.0;
            int n = NumDiffOrbs(LHS, RHS);
            if (n == 0) { H = energy(LHS); } //LHS == RHS
            else if (n == 1) //LHS(i -> a) == RHS
            {
                auto key = LHS.key() * RHS.key();
                //if element has been calculated, extract from hashtable
                if (Hia.count(key) == 1) { H = Hia.at(key); }
                else
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
                    if (si == sa) { H += pi * pa * I1(orbi, orba); } //like spin
                    //differing spin integral = 0
    
                    //two electron operator
                    std::vector<int> open, closed;
                    LHS.OpenClosed(open, closed);
                    for (int r = 0; r < closed.size(); r++)
                    {
                        int sr = closed[r] % 2;  
                        int orbr = closed[r] / 2;
    
                        if (si == sa) { H += pi * pa * I2(orbi, orba, orbr, orbr); } //direct term
                        if (si == sr && sa == sr) { H -= pi * pa * I2(orbi, orbr, orbr, orba); } //exchange term
                    }
                    //add newly calculated element to hashtable
                    Hia[key] = H;
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
                if (si == sa && sj == sb) { H += pi * pj * pa * pb * I2(orbi, orba, orbj, orbb); } //direct
                if (si == sb && sj == sa) { H -= pi * pj * pa * pb * I2(orbi, orbb, orbj, orba); } //exchange
            }
            return H;
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
            for (auto row = CI.begin(); row != CI.end(); ++row)
            {
                int j = 0;
                for (auto col = CI.begin(); col != CI.end(); ++col)
                {
                    H(i, j) = (*this)(*row, *col);
                    j++;
                }
                i++;
            }
        }

        void multiply(const FockVector &CI, const Eigen::VectorXd &z, Eigen::VectorXd &Hz) const
        {
            assert(CI.size() == z.rows());
            Hz.setZero(z.rows());
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
        }    

        /*
        void multiply(FockVector &CI) const
        {
            Eigen::VectorXd Hz = Eigen::VectorXd::Zero(CI.size());
            int i = 0;
            for (auto row = CI.begin(); row != CI.end(); ++row)
            {
                int j = 0;
                for (auto col = CI.begin(); col != CI.end(); ++col)
                {
                    Hz(i) += (*this)(*row, *col) * CI.at(*col);
                    j++;
                }
                i++;
            }
            CI.update(Hz);
            return Hz;
        }    
        */
 
    };
}
#endif
