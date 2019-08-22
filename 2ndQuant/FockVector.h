#ifndef FOCKVECTOR_HEADER
#define FOCKVECTOR_HEADER
#include "Determinant.h"
#include <iostream>
#include <unordered_set>
#include <Eigen/Dense>

struct HashDet { std::size_t operator()(const Determinant &D) const { return D.key(); } };

class FockVector
{
    protected:
    std::unordered_set<Determinant, HashDet> Store;

    public:
    std::size_t size() const { return Store.size(); } //number of determinants in vector
    auto begin() { return Store.begin(); } //begin iterator
    auto end() { return Store.end(); } //end iterator
    auto begin() const { return Store.begin(); } //begin const iterator
    auto end() const { return Store.end(); } //end const iterator
    auto insert(const Determinant &D) { return Store.emplace(D); }
    auto remove(const Determinant &D) { return Store.erase(D); }
    void clear() { Store.clear(); } //clear vector
    void ones() { for (auto it = begin(); it != end(); ++it) { it->Coeff = 1.0; } } //sets all coefficients to one
    double at(const Determinant &D) const
    {
        double val = 0.0;
        auto search = Store.find(D);
        if (search != end()) { val = search->coeff(); }
        return val;
    }
    double &at(const Determinant &D)
    {
        auto search = Store.find(D);
        if (search == end()) { search = insert(D).first; }
        return search->Coeff;
    }

    void update(const Eigen::VectorXd &V)
    {
        assert(V.size() == size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            if (std::abs(V(i)) < 1.e-8) { it->Coeff = 0.0; }
            else { it->Coeff = V(i); }
            i++;
        }
    }

    void vector(Eigen::VectorXd &V) const
    {
        V.setZero(size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            double val = it->coeff();
            if (std::abs(val) > 1.e-8) { V(i) = val; }
            i++;
        } 
    }

    void trim(double tol = 1.e-8)
    {
        auto it = begin();
        while (it != end())
        {
            if (std::abs(it->coeff()) < tol) { Store.erase(it++); }
            else { ++it; }
        } 
    }

    friend inline std::ostream &operator<<(std::ostream &os, const FockVector &V);
};
    
inline std::ostream &operator<<(std::ostream &os, const FockVector &V)
{
    for (auto it = V.begin(); it != V.end(); ++it) { std::cout << *it << std::endl; }
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
                Store.emplace(Determinant(alpha[i], beta[j]));
            }
        }
        assert(size() == alpha.size() * beta.size());
    }

    double operator()(const Determinant &D) const { return at(D); }
    double &operator()(const Determinant &D) { return at(D); }
};

class CISDVector: public FockVector
{
    private:

    public:
    CISDVector()
    {
        Determinant hf;
        hf.HartreeFock();
        std::vector<Determinant> dets;
        hf.connected(dets);
        for (int i = 0; i < dets.size(); i++) { insert(dets[i]); } 
        assert(hf.numConnected() == dets.size());
    }

    double operator()(const Determinant &D) const { return at(D); }
    double &operator()(const Determinant &D) { return at(D); }
};
#endif
