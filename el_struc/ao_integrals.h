#ifndef AO_INTEGRAL_HEADER_H
#define AO_INTEGRAL_HEADER_H
#include <iostream>
#include <Eigen/Dense>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


inline int max(int a, int b) 
{
    return a > b ? a : b;
}

inline int min(int a, int b)
{
    return a < b ? a : b;
}

class one_int
{
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version) { ar & store & norbs; }

    public:
    std::vector<double> store;
    int norbs;

    inline double &operator() (int i, int j) { return store.at(i * norbs + j); }
};

class two_int
{
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & store & zero & direct & exchange & norbs;
        }

    public:
    std::vector<double> store;
    double zero;
    int norbs;

    two_int(): zero(0.0) { }

    inline double &operator() (int i, int j, int k, int l)
    {
        zero = 0.0;
        int ij = max(i,j) * (max(i,j) + 1) / 2 + min(i,j);
        int kl = max(k,l) * (max(k,l) + 1) / 2 + min(k,l);
        int ijkl = max(ij,kl) * (max(ij,kl) + 1) / 2 + min(ij,kl);
        return store.at(ijkl);
    }
};

void read_ao_integrals(std::string fcidump, std::string s_matrix, std::string nuc, one_int &I1, two_int &I2, one_int &S, double &core_e, double &nuc_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);
#endif
