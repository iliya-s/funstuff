#ifndef INTEGRALS_HEADER_H
#define INTEGRALS_HEADER_H
#include <algorithm>
#include <vector>
#include <map>
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
        OneElectron() {}
        OneElectron(int norb) { init(norb); }
        inline void init(int norb)
        { 
            Norb = norb;
            Store.assign(Norb * Norb, 0.0);
        }

        int norb() const { return Norb; }
        inline double operator()(int i, int j) const { return Store.at(i * Norb + j); }
        inline double &operator()(int i, int j) { return Store.at(i * Norb + j); }
    };
    
    class TwoElectron
    {
        private:
        friend class boost::serialization::access;
        template <class Archive>
        void serialize(Archive &ar, const unsigned int version) { ar & Ksym & Store & Norb; }
    
        bool Ksym;
        std::vector<double> Store;
        int Norb;
    
        public:
        TwoElectron() {}
        TwoElectron(int norb, int ksym = false) { init(norb, ksym); }
        inline void init(int norb, bool ksym)
        {
            Norb = norb;
            Ksym = ksym;
            int npair = Norb * (Norb + 1) / 2;
            if (Ksym == true) { npair = Norb * Norb; }
            Store.assign(npair * (npair + 1) / 2, 0.0);
        }

        int norb() const { return Norb; }
        inline double operator()(int i, int j, int k, int l) const
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
        inline double &operator()(int i, int j, int k, int l)
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

    namespace HeatBath
    {
        class OneElectron
        {
            private:
            std::map<int, std::multimap<float, int, std::greater<float>>> Store;

            public:
            OneElectron() {}
            OneElectron(const Integral::OneElectron &I1, const Integral::TwoElectron &I2) { init(I1, I2); }
            void init(const Integral::OneElectron &I1, const Integral::TwoElectron &I2);

            auto begin(int i) { return Store.at(i).begin(); }
            auto end(int i) { return Store.at(i).end(); }
            auto begin(int i) const { return Store.at(i).begin(); }
            auto end(int i) const { return Store.at(i).end(); }            
        };
        
        class TwoElectron
        {
            private:
            std::map<std::pair<int, int>, std::multimap<float, std::pair<int, int>, std::greater<float>>> Store;

            public:
            TwoElectron() {}
            TwoElectron(const Integral::OneElectron &I1, const Integral::TwoElectron &I2) { init(I1, I2); }
            void init(const Integral::OneElectron &I1, const Integral::TwoElectron &I2);

            auto begin(std::pair<int, int> ij) { return Store.at(ij).begin(); }
            auto end(std::pair<int, int> ij) { return Store.at(ij).end(); }
            auto begin(std::pair<int, int> ij) const { return Store.at(ij).begin(); }
            auto end(std::pair<int, int> ij) const { return Store.at(ij).end(); }            
        };

    }
}

void ReadFCIDUMP(std::string FCIDUMP, Integral::OneElectron &I1, Integral::TwoElectron &I2, double &core_e, int &norb, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);

void ReadSquareMatrixIntoOneInt(std::string MatrixFile, int norb, Integral::OneElectron &I);

void ReadMatrix(std::string MatrixFile, int rows, int cols, Eigen::MatrixXd &M);

void WriteMatrix(std::string MatrixFile, const Eigen::MatrixXd &M);

void ReadAtomicOrbitalIntegrals(std::string AOFCIDUMP, std::string METRIC, Integral::OneElectron &AI1, Integral::TwoElectron &AI2, Integral::OneElectron &S, double &core_e, int &norb, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep);
#endif
