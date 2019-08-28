#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER
#include "Integrals.h"
#include "Determinant.h"
#include "FockVector.h"
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

//struct HashKey { std::size_t operator()(const std::size_t &key) const { return key; } };

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
        ar & V;
    }

    double core_e;
    int norb, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;
    Integral::OneElectron I1;
    Integral::TwoElectron I2;
    mutable std::unordered_map<std::size_t, double, HashKey> Hii;  //hashtable to store diagonal values
    mutable std::unordered_map<std::size_t, double, HashKey> Hia;  //hashtable to store singly connected values
    std::unique_ptr<FockVector> V; //current variational space of Hamiltonian

    public:
    Hamiltonian(std::string filename = "FCIDUMP")
    {
        ReadFCIDUMP(filename, I1, I2, core_e, norb, nelec, nalpha, nbeta, sz, irrep);
        InitDetVars(sz, nelec, norb);
    }

    void clear() const { Hii.clear(); Hia.clear(); } //clears hashtables 
    std::size_t size() const { return Hii.size() + Hia.size(); } //returns size of hashtables

    std::size_t dimension() const { return V->size(); } //returns dimension of variational space
    void space(const FockVector &CI) { V = std::make_unique<FockVector>(CI); }

    double energy(const Determinant &D) const
    {
        double E = 0.0;
        auto key = D.key();
        if (Hii.count(key) == 1) { E = Hii.at(key); } //if element has been calculated, extract from hashtable
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

            Hii[key] = E; //add newly calculated element to hashtable
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
            if (Hia.count(key) == 1) { H = Hia.at(key); } //if element has been calculated, extract from hashtable
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

    void diagonal(Eigen::VectorXd &CI) const
    {
        assert(V != nullptr);
        CI.setZero(V->size());
        int i = 0;
        for (auto it = V->begin(); it != V->end(); ++it)
        {
            CI(i) = energy(it->second);
            i++;
        }
    }

    void matrix(Eigen::MatrixXd &H) const
    {
        assert(V != nullptr);
        H.setZero(V->size(), V->size());
        int i = 0;
        for (auto row = V->begin(); row != V->end(); ++row)
        {
            int j = 0;
            for (auto col = V->begin(); col != V->end(); ++col)
            {
                H(i, j) = (*this)(row->second, col->second);
                j++;
            }
            i++;
        }
    }

    /*
    void multiply(const Eigen::VectorXd &z, Eigen::VectorXd &Hz) const
    {
        assert(V != nullptr);
        assert(V->size() == z.rows());
        Hz.setZero(z.rows());
        int i = 0;
        for (auto row = V->begin(); row != V->end(); ++row)
        {
            int j = 0;
            for (auto col = V->begin(); col != V->end(); ++col)
            {
                Hz(i) += (*this)(row->second, col->second) * z(j);
                j++;
            }
            i++;
        }
    }    

    void Fmultiply(const Eigen::VectorXd &z, Eigen::VectorXd &Hz) const
    {
        assert(V != nullptr);
        assert(V->size() == z.rows());
        V->update(z);
        Hz.setZero(z.rows());
        int i = 0;
        for (auto row = V->begin(); row != V->end(); ++row)
        {
            int j = 0;
            std::vector<Determinant> dets;
            row->second.connected(dets);
            for (auto col = dets.begin(); col != dets.end(); ++col)
            {
                double val = V->at(*col);
                if (val != 0.0) { Hz(i) += (*this)(row->second, *col) * val; }
                j++;
            }
            i++;
        }
        V->ones();
    }    
    */

    void multiply(const Eigen::VectorXd &z, Eigen::VectorXd &Hz) const
    {
        assert(V != nullptr);
        assert(V->size() == z.rows());
        Hz.setZero(z.rows());
        Determinant det;
        det.HartreeFock();
    
        if (V->size() <= det.numConnected()) //if the fock space is on the order of the number of connections, normal matrix multiplication
        {
            int i = 0;
            for (auto row = V->begin(); row != V->end(); ++row)
            {
                int j = 0;
                for (auto col = V->begin(); col != V->end(); ++col)
                {
                    Hz(i) += (*this)(row->second, col->second) * z(j);
                    j++;
                }
                i++;
            }
        }    
        else //if the fock space is larger than the number of connections, fast matrix multiplication
        {
            V->update(z);
            int i = 0;
            for (auto row = V->begin(); row != V->end(); ++row)
            {
                std::vector<Determinant> dets;
                row->second.connected(dets);
                for (auto col = dets.begin(); col != dets.end(); ++col)
                {
                    double val = V->at(*col);
                    if (val != 0.0) { Hz(i) += (*this)(row->second, *col) * val; }
                }
                i++;
            }
            V->ones();
        }
    }

    Eigen::VectorXd operator*(const Eigen::VectorXd &z) const
    {
        assert(V != nullptr);
        assert(V->size() == z.rows());
        Eigen::VectorXd Hz = Eigen::VectorXd::Zero(z.rows());
        Determinant det;
        det.HartreeFock();
    
        if (V->size() <= det.numConnected()) //if the fock space is on the order of the number of connections, normal matrix multiplication
        {
            int i = 0;
            for (auto row = V->begin(); row != V->end(); ++row)
            {
                int j = 0;
                for (auto col = V->begin(); col != V->end(); ++col)
                {
                    Hz(i) += (*this)(row->second, col->second) * z(j);
                    j++;
                }
                i++;
            }
        }    
        else //if the fock space is larger than the number of connections, fast matrix multiplication
        {
            V->update(z);
            int i = 0;
            for (auto row = V->begin(); row != V->end(); ++row)
            {
                std::vector<Determinant> dets;
                row->second.connected(dets);
                for (auto col = dets.begin(); col != dets.end(); ++col)
                {
                    double val = V->at(*col);
                    if (val != 0.0) { Hz(i) += (*this)(row->second, *col) * val; }
                }
                i++;
            }
            V->ones();
        }
        return Hz;
    }


};

#endif
