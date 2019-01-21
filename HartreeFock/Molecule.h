#ifndef MOLECULE_HEADER_H
#define MOLECULE_HEADER_H
#include "AtomicIntegrals.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
 
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
            ar & nuc_e;
            ar & norbs;
            ar & nelec;
            ar & nalpha;
            ar & nbeta;
            ar & sz;
            ar & irrep
            ar & AI1;
            ar & AI2;
            ar & S;
        }

    public:
    //molcule information
    double core_e, nuc_e;
    int norbs, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;
    bool built = false;
    
    //ao integrals
    AOneInt AI1, S;
    ATwoInt AI2;

    //member functions
    int build(std::string AOfcidump, std::string metric_file, std::string nuc_file)
    {
        built = true;
        return ReadAOIntegrals(AOfcidump, metric_file, nuc_file, AI1, AI2, S, core_e, nuc_e, norbs, nelec, nalpha, nbeta, sz, irrep);
    }

    //constructors
    Molecule() {}
    Molecule(std::string AOfcidump, std::string metric_file, std::string nuc_file)
    {
        build(AOfcidump, metric_file, nuc_file);
    }
}; //Molecule
#endif
