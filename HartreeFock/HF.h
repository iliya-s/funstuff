#ifndef HF_HEADER_H
#define HF_HEADER_H
#include "Molecule.h"
#include "AtomicIntegrals.h"
#include <Eigen/Dense>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

class RHF
{
    private:
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & E;
            ar & orbitals;
            ar & orbitals_e;
        }

    public:
    //hartreefock energy
    double E;
    //molecular orbitals
    Eigen::MatrixXd orbitals;
    Eigen::VectorXd orbitals_e;

    //constructors
    RHF() {}
    RHF(Molecule &mol);

    //member functions
    int Run(Molecule &mol);
}; //RHF

class UHF
{
    private:
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar & E;
            ar & orbitals;
            ar & orbitals_e;
        }

    public:
    //hartreefock energy
    double E;
    //molecular orbitals, 0 indexes alpha electrons, 1 indexes beta electrons
    std::array<Eigen::MatrixXd, 2> orbitals;
    std::array<Eigen::VectorXd, 2> orbitals_e;

    //constructors
    UHF() {}
    UHF(Molecule &mol);

    //member functions
    int Run(Molecule &mol);
}; //UHF
#endif
