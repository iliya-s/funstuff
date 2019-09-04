#ifndef FOCKVECTOR_HEADER
#define FOCKVECTOR_HEADER
#include "Determinant.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

class FockVector
{
    protected:
    std::vector<Determinant> Store;

    public:
    FockVector() {}
    std::size_t size() const { return Store.size(); } //number of determinants in vector

    auto begin() { return Store.begin(); } //begin iterator
    auto end() { return Store.end(); } //end iterator
    auto begin() const { return Store.begin(); } //begin const iterator
    auto end() const { return Store.end(); } //end const iterator

    void insert(const Determinant &D)
    {
        auto it = std::lower_bound(begin(), end(), D);
        if (D == *it) { return; }
        else { Store.insert(it, D); }
    }
    void clear() { Store.clear(); } //clear vector
    void ones() { for (auto it = begin(); it != end(); ++it) { it->coeff(1.0); } } //sets all coefficients to one

    int at(const Determinant &D) const
    {
        int i = size();
        auto it = std::lower_bound(begin(), end(), D);
        if (D == *it) { i = std::distance(begin(), it); }
        return i;
    }
    int operator()(const Determinant &D) const { return at(D); }

    const Determinant &at(int i) const { return Store.at(i); }
    const Determinant &operator()(int i) const { return at(i); }

    void update(const Eigen::VectorXd &V)
    {
        assert(V.size() == size());
        for (int i = 0; i < size(); i++) { Store.at(i).coeff(V(i)); } 
    }

    void vector(Eigen::VectorXd &V) const
    {
        V.setZero(size());
        for (int i = 0; i < size(); i++) { V(i) = Store.at(i).coeff(); } 
    }

    void determinants(std::vector<Determinant> &dets)
    {
        dets.clear();
        for (auto it = begin(); it != end(); ++it) { dets.push_back(*it); }
        std::sort(dets.begin(), dets.end(), [](const Determinant &a, const Determinant &b) { return  a.coeff() > b.coeff(); });
    }

    void trim(double tol = 1.e-4)
    {
        for(auto it = begin(); it != end();)
        {
            if(std::abs(it->coeff()) < tol) { Store.erase(it); }
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
                Determinant det(alpha[i], beta[j]);
                Store.push_back(det);
            }
        }
        std::sort(begin(), end());
        assert(Store.size() == alpha.size() * beta.size());
    }
};

class CISVector: public FockVector
{
    private:

    public:
    CISVector()
    {
        Determinant hf;
        hf.HartreeFock();
        Store.push_back(hf);
        std::vector<Determinant> dets;
        hf.singleExcitations(dets);
        for (int i = 0; i < dets.size(); i++)
        {
            Store.push_back(dets[i]);
        }
        std::sort(begin(), end());
        assert(hf.nSingleExcitations() + 1 == Store.size());
    }
};

class CIDVector: public FockVector
{
    private:

    public:
    CIDVector()
    {
        Determinant hf;
        hf.HartreeFock();
        Store.push_back(hf);
        std::vector<Determinant> dets;
        hf.doubleExcitations(dets);
        for (int i = 0; i < dets.size(); i++)
        {
            Store.push_back(dets[i]);
            std::sort(begin(), end());
        }
        std::sort(begin(), end());
        assert(hf.nDoubleExcitations() + 1 == Store.size());
    }
};

class CISDVector: public FockVector
{
    private:

    public:
    CISDVector()
    {
        Determinant hf;
        hf.HartreeFock();
        Store.push_back(hf);
        std::vector<Determinant> dets;
        hf.excitations(dets);
        for (int i = 0; i < dets.size(); i++)
        {
            Store.push_back(dets[i]);
            std::sort(begin(), end());
        }
        std::sort(begin(), end());
        assert(hf.nExcitations() + 1 == Store.size());
    }
};
#endif
