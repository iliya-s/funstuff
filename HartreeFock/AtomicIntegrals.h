#ifndef AO_INTEGRAL_HEADER_H
#define AO_INTEGRAL_HEADER_H
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

class AOneInt
{
    private:
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & store & norbs;
        }

    public:
    std::vector<double> store;
    int norbs;

    inline double operator() (int i, int j) const
    {
        return store.at(i * norbs + j);
    }
    inline double &operator() (int i, int j)
    {
        return store.at(i * norbs + j);
    }
};

class ATwoInt
{
    private:
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & store & norbs;
        }

    public:
    std::vector<double> store;
    int norbs;

    inline double operator() (int i, int j, int k, int l) const
    {
        int ij = std::max(i, j) * (std::max(i, j) + 1) / 2 + std::min(i, j);
        int kl = std::max(k, l) * (std::max(k, l) + 1) / 2 + std::min(k, l);
        int ijkl = std::max(ij, kl) * (std::max(ij, kl) + 1) / 2 + std::min(ij, kl);
        return store.at(ijkl);
    }
    inline double &operator() (int i, int j, int k, int l)
    {
        int ij = std::max(i, j) * (std::max(i, j) + 1) / 2 + std::min(i, j);
        int kl = std::max(k, l) * (std::max(k, l) + 1) / 2 + std::min(k, l);
        int ijkl = std::max(ij, kl) * (std::max(ij, kl) + 1) / 2 + std::min(ij, kl);
        return store.at(ijkl);
    }
};

int ReadAOIntegrals(std::string AOfcidump, std::string metric, std::string nuc, AOneInt &AI1, ATwoInt &AI2, AOneInt &S, double &core_e, double &nuc_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);
#endif
