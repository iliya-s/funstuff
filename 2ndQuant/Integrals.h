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
        void serialize(Archive &ar, const unsigned int version) { ar & Store & Norb; }
    
        std::vector<double> Store;
        int Norb;
    
        public:
        inline void init(int norb)
        { 
            Norb = norb;
            Store.resize(Norb * Norb);
        }

        inline double operator() (int i, int j) const { return Store.at(i * Norb + j); }
        inline double &operator() (int i, int j) { return Store.at(i * Norb + j); }
    };
    
    class TwoElectron
    {
        private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version) { ar & Store & Norb; }
    
        bool Ksym;
        std::vector<double> Store;
        int Norb;
    
        public:
        inline void init(int norb, bool ksym)
        {
            Norb = norb;
            Ksym = ksym;
            int npair = Norb * (Norb + 1) / 2;
            if (Ksym == true)
            {
                npair = Norb * Norb;
            }
            Store.resize(npair * (npair + 1) / 2);
        }

        inline double operator() (int i, int j, int k, int l) const
        {
            int ij, kl;
            if (Ksym == true)
            {
                ij = i * Norb + j;
                kl = k * Norb + l;
            }
            else
            {
                ij = std::max(i, j) * (std::max(i, j) + 1) / 2 + std::min(i, j);
                kl = std::max(k, l) * (std::max(k, l) + 1) / 2 + std::min(k, l);
            }
            int ijkl = std::max(ij, kl) * (std::max(ij, kl) + 1) / 2 + std::min(ij, kl);
            return Store.at(ijkl);
        }
        inline double &operator() (int i, int j, int k, int l)
        {
            int ij, kl;
            if (Ksym == true)
            {
                ij = i * Norb + j;
                kl = k * Norb + l;
            }
            else
            {
                ij = std::max(i, j) * (std::max(i, j) + 1) / 2 + std::min(i, j);
                kl = std::max(k, l) * (std::max(k, l) + 1) / 2 + std::min(k, l);
            }
            int ijkl = std::max(ij, kl) * (std::max(ij, kl) + 1) / 2 + std::min(ij, kl);
            return Store.at(ijkl);
        }
    };
}

void ReadFCIDUMP(std::string FCIDUMP, Integral::OneElectron &I1, Integral::TwoElectron &I2, double &core_e, int &Norb, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);

void ReadSquareMatrixIntoOneInt(std::string MatrixFile, int Norb, Integral::OneElectron &I);

void ReadMatrix(std::string MatrixFile, int rows, int cols, Eigen::MatrixXd &M);

void WriteMatrix(std::string MatrixFile, const Eigen::MatrixXd &M);

void ReadAtomicOrbitalIntegrals(std::string AOFCIDUMP, std::string METRIC, Integral::OneElectron &AI1, Integral::TwoElectron &AI2, Integral::OneElectron &S, double &core_e, int &Norb, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);
#endif
