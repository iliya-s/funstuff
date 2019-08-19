#ifndef FOCKVECTOR_HEADER
#define FOCKVECTOR_HEADER
#include "Determinant.h"
#include <map>
#include <iostream>
#include <Eigen/Dense>

class FockVector
{
    protected:
    std::map<std::size_t, Determinant> Store;

    public:
    std::size_t size() const { return Store.size(); }
    auto begin() { return Store.begin(); }
    auto end() { return Store.end(); }
    auto begin() const { return Store.begin(); }
    auto end() const { return Store.end(); }
    double at(const Determinant &D) const { return Store.at(D.key()).coeff(); }
    double &at(const Determinant &D) { return Store.at(D.key()).coeff(); }

    void update(const Eigen::VectorXd &V)
    {
        assert(V.size() == size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            if (std::abs(V(i)) < 1.e-8) { it->second.coeff(0.0); }
            else { it->second.coeff(V(i)); }
            i++;
        }
    }

    void vector(Eigen::VectorXd &V) const
    {
        V.setZero(size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            double val = it->second.coeff();
            if (std::abs(val) > 1.e-8) { V(i) = val; }
            i++;
        } 
    }

    void trim(double e = 1.e-8)
    {
        auto it = begin();
        while (it != end())
        {
            if (std::abs(it->second.coeff()) < e) { Store.erase(it++); }
            else { ++it; }
        } 
    }

    friend inline std::ostream &operator<<(std::ostream &os, const FockVector &V);
};
    
inline std::ostream &operator<<(std::ostream &os, const FockVector &V)
{
    for (auto it = V.begin(); it != V.end(); ++it)
    {
        std::cout << it->second << std::endl;
    }
    return os;
}

class FCIVector: public FockVector
{
    private:

    public:
    FCIVector()
    {
        std::vector<std::vector<int>> alpha, beta;
        GenerateCombinations(Determinant::norb(), Determinant::nalpha(), alpha);
        GenerateCombinations(Determinant::norb(), Determinant::nbeta(), beta);
        for (int i = 0; i < alpha.size(); i++)
        {
            for (int j = 0; j < beta.size(); j++)
            {
                Determinant det(alpha[i], beta[j]);
                Store.emplace(det.key(), det);
            }
        }
        assert(size() == alpha.size() * beta.size());
    }

    double operator()(const Determinant &D) const { return at(D); }
    double &operator()(const Determinant &D) { return at(D); }
};
#endif
