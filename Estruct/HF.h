#ifndef HF_HEADER_H
#define HF_HEADER_H
#include "Molecule.h"
#include "Integral.h"
#include <Eigen/Dense>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace boost;
using namespace Eigen;

class RHF
{
    friend class serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar
            & E
            & orbitals
            & orbitals_e;
        }

    public:
    //hartreefock energy
    double E;
    //molecular orbitals
    MatrixXd orbitals;
    VectorXd orbitals_e;

    //constructors
    RHF() {}
    RHF(Molecule &mol);

    //member functions
    int run(Molecule &mol);
};

class UHF
{
    friend class serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar
            & E
            & alpha_orbitals
            & beta_orbitals
            & alpha_orbitals_e
            & beta_orbitals_e;
        }

    public:
    //hartreefock energy
    double E;
    //molecular orbitals
    MatrixXd alpha_orbitals;
    MatrixXd beta_orbitals;
    VectorXd alpha_orbitals_e;
    VectorXd beta_orbitals_e;

    //constructors
    UHF() {}
    UHF(Molecule &mol);

    //member functions
    int run(Molecule &mol);
};
#endif
