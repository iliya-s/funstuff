#ifndef MOLECULE_HEADER_H
#define MOLECULE_HEADER_H
#include "Integral.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;
using namespace boost;

class Molecule
{
    friend class serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar
            & built
            & core_e
            & nuc_e
            & norbs
            & nelec
            & nalpha
            & nbeta
            & sz
            & irrep
            & AI1
            & AI2
            & M;
        }

    public:
    //molcule information
    double core_e, nuc_e;
    int norbs, nelec, nalpha, nbeta, sz;
    vector<int> irrep;
    bool built = false;
    
    //ao integrals
    OneInt AI1, M;
    TwoInt AI2;

    //constructors
    Molecule() {}
    Molecule(string ao_fcidump, string metric_file, string nuc_file);

    //member functions
    int build(string ao_fcidump, string metric_file, string nuc_file);

    /*
    void integral_transform(string hf);
    */
};
#endif
