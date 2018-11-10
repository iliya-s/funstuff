#ifndef INTEGRAL_HEADER_H
#define INTEGRAL_HEADER_H
#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace std;
using namespace boost;
using namespace Eigen;

class OneInt
{
    friend class serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & store & norbs;
        }

    public:
    MatrixXd store;
    int norbs;

    inline double &operator() (int i, int j)
    {
        return store(i,j);
    }
};

class TwoInt
{
    friend class serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & store & norbs;
        }

    public:
    vector<double> store;
    int norbs;

    inline double &operator() (int i, int j, int k, int l)
    {
        int ij = max(i,j) * (max(i,j) + 1) / 2 + min(i,j);
        int kl = max(k,l) * (max(k,l) + 1) / 2 + min(k,l);
        int ijkl = max(ij,kl) * (max(ij,kl) + 1) / 2 + min(ij,kl);
        return store.at(ijkl);
    }
};

int read_ao_integrals(string fcidump, string s_matrix, string nuc, OneInt &AI1, TwoInt &AI2, OneInt &S, double &core_e, double &nuc_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, vector<int> &irrep);
#endif
