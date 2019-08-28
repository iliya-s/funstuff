#ifndef FOCKVECTOR_HEADER
#define FOCKVECTOR_HEADER
#include "Determinant.h"
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

struct HashKey { std::size_t operator()(const std::size_t &key) const { return key; } };

class FockVector
{
    protected:
    std::unordered_map<std::size_t, Determinant, HashKey> Store;

    public:
    FockVector() {}
    std::size_t size() const { return Store.size(); } //number of determinants in vector
    auto begin() { return Store.begin(); } //begin iterator
    auto end() { return Store.end(); } //end iterator
    auto begin() const { return Store.begin(); } //begin const iterator
    auto end() const { return Store.end(); } //end const iterator
    void insert(const Determinant &D)
    {
        auto pair = Store.emplace(D.key(), D);
        if (pair.second == false)
        {
            std::cout << "failed insertion" << std::endl;
            std::cout << D.key() << " " << D << std::endl;
            std::cout << pair.first->second.key() << " " << pair.first->second << std::endl;
        }
    }
    auto remove(const Determinant &D) { return Store.erase(D.key()); }
    void clear() { Store.clear(); } //clear vector
    void ones() { for (auto it = begin(); it != end(); ++it) { it->second.coeff(1.0); } } //sets all coefficients to one
    double at(const Determinant &D) const
    {
        double val = 0.0;
        if (Store.count(D.key()) == 1) { val = Store.at(D.key()).coeff(); }
        return val;
    }

    double operator()(const Determinant &D) const { return at(D); }

    void update(const Eigen::VectorXd &V)
    {
        assert(V.size() == size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            it->second.coeff(V(i));
            i++;
        }
    }

    void vector(Eigen::VectorXd &V) const
    {
        V.setZero(size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            V(i) = it->second.coeff();
            i++;
        } 
    }

    void determinants(std::vector<Determinant> &dets)
    {
        dets.clear();
        for (auto it = begin(); it != end(); ++it)
        {
            dets.push_back(it->second);
        }
        std::sort(dets.begin(), dets.end(), [](const Determinant &a, const Determinant &b) { return  a.coeff() > b.coeff(); });
    }

    void trim(double tol = 1.e-4)
    {
        auto it = begin();
        while (it != end())
        {
            if (std::abs(it->second.coeff()) < tol) { Store.erase(it++); }
            else { ++it; }
        } 
    }

    friend inline std::ostream &operator<<(std::ostream &os, const FockVector &V);
};
    
inline std::ostream &operator<<(std::ostream &os, const FockVector &V)
{
    for (auto it = V.begin(); it != V.end(); ++it) { std::cout << it->second << std::endl; }
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
                insert(det);
            }
        }
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
        insert(hf);
        std::vector<Determinant> dets;
        hf.singlyConnected(dets);
        for (int i = 0; i < dets.size(); i++) { insert(dets[i]); } 
        assert((hf.numSinglyConnected() + 1) == Store.size());
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
        insert(hf);
        std::vector<Determinant> dets;
        hf.doublyConnected(dets);
        for (int i = 0; i < dets.size(); i++) { insert(dets[i]); } 
        assert((hf.numDoublyConnected() + 1) == Store.size());
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
        std::vector<Determinant> dets;
        hf.connected(dets);
        for (int i = 0; i < dets.size(); i++) { insert(dets[i]); } 
        assert(hf.numConnected() == Store.size());
    }
};
#endif
