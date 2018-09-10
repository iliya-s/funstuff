#ifndef MOLECULE_HEADER_H
#define MOLECULE_HEADER_H
#include "integrals.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
 
#ifndef SERIAL
#include <boost/mpi.hpp>
#endif


class molecule
{
    friend class boost::serialization::access;
    template <class Archive>
        void serialize(Archive &ar, const unsigned int version)
        {
            ar
            & core_e
            & built
            & nuc_e
            & norbs
            & nelec
            & nalpha
            & nbeta
            & sz
            & irrep
            & ao_i1
            & ao_i2
            & metric
            & E
            & mo_i1
            & mo_i2
            & orbitals
            & orbitals_e
            & alpha_orbitals
            & alpha_orbitals_e
            & beta_orbitals
            & beta_orbitals_e;
        }

    public:
    //molcule information
    double core_e, nuc_e;
    int norbs, nelec, nalpha, nbeta, sz;
    std::vector<int> irrep;
    bool built = false;
    
    //ao integrals
    one_int ao_i1, metric;
    two_int ao_i2;

    //hartree fock variables
    double E;
    one_int mo_i1;
    two_int mo_i2;
        //used for restricted hartree fock
    Eigen::MatrixXd orbitals;
    Eigen::VectorXd orbitals_e;
        //used for unrestriced hartree fock
    Eigen::MatrixXd alpha_orbitals;
    Eigen::VectorXd alpha_orbitals_e;
    Eigen::MatrixXd beta_orbitals;
    Eigen::VectorXd beta_orbitals_e;

    //member functions
    void build(std::string ao_integrals, std::string metric_file, std::string nuc);
    void rhf();
    void uhf();
};

#endif
