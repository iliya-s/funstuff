#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER
#include "Integrals.h"
#include "Determinant.h"
#include "FockVector.h"
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

struct HashKey { std::size_t operator()(const std::size_t &key) const { return key; } };

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
        ar & HBI1 & HBI2;
        ar & Hii & Hia;
        ar & V;
    }

    double core_e;
    int norb, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;
    Integral::OneElectron I1;
    Integral::TwoElectron I2;
    Integral::HeatBath::OneElectron HBI1;
    Integral::HeatBath::TwoElectron HBI2;
    mutable std::unordered_map<std::size_t, double, HashKey> Hii;  //hashtable to store diagonal values
    mutable std::unordered_map<std::size_t, double, HashKey> Hia;  //hashtable to store singly connected values
    std::unique_ptr<FockVector> V; //current variational space of Hamiltonian

    public:
    //constructor
    Hamiltonian(std::string filename = "FCIDUMP")
    {
        ReadFCIDUMP(filename, I1, I2, core_e, norb, nelec, nalpha, nbeta, sz, irrep);
        InitDetVars(sz, nelec, norb);
        HBI1.init(I1, I2);
        HBI2.init(I1, I2);
    }

    //returns integrals
    const Integral::OneElectron &oneElectron() const { return I1; }
    const Integral::TwoElectron &twoElectron() const { return I2; }
    const Integral::HeatBath::OneElectron &heatBathOneElectron() const { return HBI1; }
    const Integral::HeatBath::TwoElectron &heatBathTwoElectron() const { return HBI2; }

    //functions for Hamiltonian hashtable
    void clear() const { Hii.clear(); Hia.clear(); } //clears hashtables 
    std::size_t size() const { return Hii.size() + Hia.size(); } //returns size of hashtables

    //functions for variational space
    std::size_t dimension() const { return V->size(); } //returns dimension of variational space
    void space(const FockVector &CI) { V = std::make_unique<FockVector>(CI); } //sets variational space

    //energy of determinant
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

    //matrix element between two determinants
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
                auto ia = OneDiffOrbIndices(LHS, RHS);
                int i = ia.first;
                int a = ia.second;

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

                    if (si == sa) { H += pi * pa * I2(orbi, orba, orbr, orbr); }
                    if (si == sr && sa == sr) { H -= pi * pa * I2(orbi, orbr, orbr, orba); }
                }
                //add newly calculated element to hashtable
                Hia[key] = H;
            }
        }
        else if (n == 2) //LHS(ij -> ab) == RHS
        {
            auto ijab = TwoDiffOrbIndices(LHS, RHS);
            int i = ijab.first.first;
            int j = ijab.first.second;
            int a = ijab.second.first;
            int b = ijab.second.second;

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
            if (si == sa && sj == sb) { H += pi * pj * pa * pb * I2(orbi, orba, orbj, orbb); }
            if (si == sb && sj == sa) { H -= pi * pj * pa * pb * I2(orbi, orbb, orbj, orba); }
        }
        return H;
    }

    void diagonal(Eigen::VectorXd &CI) const
    {
        assert(V != nullptr);
        CI.setZero(V->size());
        for (int i = 0; i < CI.size(); i++) { CI(i) = energy(V->at(i)); }
    }

    void matrix(Eigen::MatrixXd &H) const
    {
        assert(V != nullptr);
        H.setZero(V->size(), V->size());
        for (int i = 0; i < V->size(); i++)
        {
            for (int j = 0; j < V->size(); j++)
            {
                H(i, j) = (*this)(V->at(i), V->at(j));
            }
        }
    }

    void multiply(const Eigen::VectorXd &z, Eigen::VectorXd &Hz) const
    {
        assert(V != nullptr);
        assert(V->size() == z.rows());
        Hz.setZero(z.rows());
        Determinant det;
        det.HartreeFock();

        if (V->size() <= det.nExcitations() + 1) //if the fock space is on the order of the number of connections, normal matrix multiplication
        //if (true)
        {
            for (int i = 0; i < V->size(); i++)
            {
                for (int j = 0; j < V->size(); j++)
                {
                    Hz(i) += (*this)(V->at(i), V->at(j)) * z(j);
                }
            }
        }    
        else //if the fock space is larger than the number of connections, fast matrix multiplication
        {
            for (int i = 0; i < V->size(); i++)
            {
                const Determinant &det(V->at(i));
                //std::vector<Determinant> dets;
                //dets.push_back(det);
                //det.excitations(dets);
                //det.screenedExcitations(HBI1, HBI2, dets);
                //for (int j = 0; j < dets.size(); j++)
                //{
                //    int nj = V->at(dets[j]);
                //    if (nj < V->size()) { Hz(i) += (*this)(det, dets[j]) * z(nj); }
                //}
                std::vector<std::pair<int, int>> s;
                det.singleExcitations(s);
                std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> d;
                det.doubleExcitations(d);
                //self
                Hz(i) += energy(det) * z(i);
                //single excitations
                for (int j = 0; j < s.size(); j++)
                {
                    int &p = s[j].first;
                    int &q = s[j].second;
                    Determinant dexcite(det);
                    dexcite.set(p, false);
                    dexcite.set(q, true);
                    int nexcite = V->at(dexcite);
                    if (nexcite < V->size()) { Hz(i) += (*this)(det, dexcite) * z(nexcite); }
                }
                //double excitations
                for (int j = 0; j < d.size(); j++)
                {
                    int &p = d[j].first.first;
                    int &q = d[j].first.second;
                    int &r = d[j].second.first;
                    int &s = d[j].second.second;
                    Determinant dexcite(det);
                    dexcite.set(p, false);
                    dexcite.set(q, false);
                    dexcite.set(r, true);
                    dexcite.set(s, true);
                    int nexcite = V->at(dexcite);
                    if (nexcite < V->size()) { Hz(i) += (*this)(det, dexcite) * z(nexcite); }
                }
            }
        }
    }

    Eigen::VectorXd operator*(const Eigen::VectorXd &z) const
    {
        assert(V != nullptr);
        assert(V->size() == z.rows());
        Eigen::VectorXd Hz = Eigen::VectorXd::Zero(z.rows());
        Determinant det;
        det.HartreeFock();

        if (V->size() <= det.nExcitations() + 1) //if the fock space is on the order of the number of connections, normal matrix multiplication
        {
            for (int i = 0; i < V->size(); i++)
            {
                for (int j = 0; j < V->size(); j++)
                {
                    Hz(i) += (*this)(V->at(i), V->at(j)) * z(j);
                }
            }
        }    
        else //if the fock space is larger than the number of connections, fast matrix multiplication
        {
            for (int i = 0; i < V->size(); i++)
            {
                const Determinant &det(V->at(i));
                std::vector<Determinant> dets;
                dets.push_back(det);
                det.screenedExcitations(HBI1, HBI2, dets, 0.0);
                for (int j = 0; j < dets.size(); j++)
                {
                    int nj = V->at(dets[j]);
                    if (nj < V->size()) { Hz(i) += (*this)(det, dets[j]) * z(nj); }
                }
            }
        }
        return Hz;
    }


};
#endif
