#ifndef INTEGRALS_HEADER_H
#define INTEGRALS_HEADER_H
#include <algorithm>
#include <Eigen/Dense>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace Integral
{
    class OneElectron
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
    
    class TwoElectron
    {
        private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & store & norbs;
        }
    
        public:
        bool ksym;
        std::vector<double> store;
        int norbs;
        Eigen::MatrixXd Direct, Exchange;
    
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
}

void ReadFCIDUMP(std::string FCIDUMP, Integral::OneElectron &I1, Integral::TwoElectron &I2, double &core_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);

void ReadSquareMatrixIntoOneInt(std::string MatrixFile, int norbs, Integral::OneElectron &I);

void ReadMatrix(std::string MatrixFile, int rows, int cols, Eigen::MatrixXd &M);

void WriteMatrix(std::string MatrixFile, const Eigen::MatrixXd &M);

void ReadAtomicOrbitalIntegrals(std::string AOFCIDUMP, std::string METRIC, Integral::OneElectron &AI1, Integral::TwoElectron &AI2, Integral::OneElectron &S, double &core_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);
#endif
