#include "Determinant.h"
#include <map>
#include <iostream>

class CIVector
{
    protected:
    std::map<std::size_t, Determinant> Store;

    public:
    std::size_t size() { return Store.size(); }
    auto begin() const { return Store.begin(); }
    auto end() const { return Store.end(); }
    double at(const Determinant &D) const { return Store.at(D.key()).coeff(); }
    double &at(const Determinant &D) { return Store.at(D.key()).coeff(); }
    friend std::ostream &operator<<(std::ostream &os, const CIVector &V);
};
    
std::ostream &operator<<(std::ostream &os, const CIVector &V)
{
    for (auto it = V.begin(); it != V.end(); ++it)
    {
        std::cout << it->second << std::endl;
    }
    return os;
}

class FCIVector: public CIVector
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
        assert(Store.size() == alpha.size() * beta.size());
    }

    Determinant operator()(int i) //this is kind of expensive, std::advance() is linear in i
    {
        assert(i >= 0 && i < Store.size());
        auto it = begin();
        std::advance(it, i);
        return it->second;
    }
    
    double operator()(const Determinant &D) const
    {
        return at(D);
    }
    double &operator()(const Determinant &D)
    {
        return at(D);
    }
};
