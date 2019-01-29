#ifndef MOLECULE_HEADER_H
#define MOLECULE_HEADER_H
#include "Integrals.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

class Molecule
{
    private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & built;
        ar & core_e;
        ar & norbs;
        ar & nelec;
        ar & nalpha;
        ar & nbeta;
        ar & sz;
        ar & irrep
        ar & I1;
        ar & I2;
        ar & S;
    }

    public:
    //molcule information
    bool built = false;
    double core_e;
    int norbs, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;
    
    //atomic orbital integrals
    Integral::OneElectron I1, S;
    Integral::TwoElectron I2;

    //member functions
    void build(std::string AOFCIDUMP, std::string METRIC)
    {
        built = true;
        ReadAtomicOrbitalIntegrals(AOFCIDUMP, METRIC, I1, I2, S, core_e, norbs, nelec, nalpha, nbeta, sz, irrep);
    }

    //constructors
    Molecule() {}
    Molecule(std::string AOFCIDUMP, std::string METRIC)
    {
        build(AOFCIDUMP, METRIC);
    }
}; //Molecule
#endif
