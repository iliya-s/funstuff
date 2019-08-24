#ifndef FOCKVECTOR_HEADER
#define FOCKVECTOR_HEADER
#include "Determinant.h"
#include <iostream>
#include <unordered_map>
#include <Eigen/Dense>

struct HashDet { std::size_t operator()(const Determinant &D) const { return D.key(); } };
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
    auto insert(const Determinant &D) { return Store.emplace(D.key(), D); }
    auto remove(const Determinant &D) { return Store.erase(D.key()); }
    void clear() { Store.clear(); } //clear vector
    void ones() { for (auto it = begin(); it != end(); ++it) { it->second.coeff(1.0); } } //sets all coefficients to one
    double at(const Determinant &D) const
    {
        double val = 0.0;
        //auto search = Store.find(D);
        if (Store.count(D.key()) == 1) { val = Store.at(D.key()).coeff(); }

        return val;
    }
    /*
    double &at(const Determinant &D)
    {
        std::cout << "i" << std::endl;
        auto search = Store.find(D);
        if (search == end()) { search = insert(D).first; }
        return search->Coeff;
    }
    */
    double operator()(const Determinant &D) const { return at(D); }
    //double &operator()(const Determinant &D) { return at(D); }

    void update(const Eigen::VectorXd &V)
    {
        assert(V.size() == size());
        int i = 0;
        for (auto it = begin(); it != end(); ++it)
        {
            //if (std::abs(V(i)) < 1.e-8) { it->Coeff = 0.0; }
            //else { it->Coeff = V(i); }
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
            //double val = it->coeff();
            //if (std::abs(val) > 1.e-8) { V(i) = val; }
            V(i) = it->second.coeff();
            i++;
        } 
    }

    void trim(double tol = 1.e-8)
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
                Store.emplace(det.key(), det);
            }
        }
        assert(size() == alpha.size() * beta.size());
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
        assert(hf.numConnected() == dets.size());
    }
};
#endif
