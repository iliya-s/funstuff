#include "ao_integrals.h"
#include <fstream>
#include <math.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string.hpp>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

void read_ao_integrals(std::string fcidump, std::string s_matrix, std::string nuc, one_int &I1, two_int &I2, one_int &S, double &core_e, double &nuc_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep)
{
    std::ifstream f1(fcidump);
    if (!f1.is_open())
    {
        std::cout << "Integral file " << fcidump << " does not exist." << std::endl;
        exit(1);
    }
    std::ifstream f2(s_matrix);
    if (!f2.is_open())
    {
        std::cout << "Integral file " << s_matrix << " does not exist." << std::endl;
        exit(2);
    }
    std::ifstream f3(nuc);
    if (!f3.is_open())
    {
        std::cout << "Nuclear energy file " << nuc << " does not exist." << std::endl;
    }

    bool start_scaling = false;
    norbs = -1;
    nelec = -1;
    sz = -1;

    int index = 0;
    std::vector<std::string> tok;
    std::string msg;
    while (!f1.eof())
    {
        std::getline(f1, msg);
        boost::algorithm::trim(msg);
        boost::split(tok, msg, boost::is_any_of(", \t="), boost::token_compress_on);
        if (start_scaling == false && tok.size() == 1 && (boost::iequals(tok[0], "&END") || boost::iequals(tok[0], "/")))
        {
            start_scaling = true;
            index += 1;
            break;
        }
        else if (start_scaling == false)
        {
            if (boost::iequals(tok[0].substr(0,4), "&FCI"))
            {
                if (boost::iequals(tok[1].substr(0,4), "NORB")) { norbs = atoi(tok[2].c_str()); }
                if (boost::iequals(tok[3].substr(0,5), "NELEC")) { nelec = atoi(tok[4].c_str()); }
                if (boost::iequals(tok[5].substr(0,5), "MS2")) { sz = atoi(tok[6].c_str()); }
            } 
            else if (boost::iequals(tok[0].substr(0,4), "ISYM")) { continue; }
            else if (boost::iequals(tok[0].substr(0,4), "KSYM")) { continue; }
            else if (boost::iequals(tok[0].substr(0,6), "ORBSYM"))
            {
                for (int i = 1; i < tok.size(); i++) { irrep.push_back(atoi(tok[i].c_str())); }
            }
            else 
            {
                for (int i = 0; i < tok.size(); i++) { irrep.push_back(atoi(tok[i].c_str())); }
            }
            index += 1;
        }
    } //while

    if (norbs == -1 || nelec == -1 || sz == -1)
    {
        std::cout << "failed reading norbs, nelec, or ms2" << std::endl;
        exit(3);
    }
    nalpha = nelec / 2 + sz;
    nbeta = nelec - nalpha;
    irrep.resize(norbs);
    core_e = 0.0;

    I2.store.clear();
    I2.norbs = norbs;
    I2.store.resize((norbs * norbs * norbs * norbs));

    I1.store.clear();
    I1.norbs = norbs;
    I1.store.resize(norbs * norbs, 0.0);

    S.store.clear();
    S.norbs = norbs;
    S.store.resize(norbs * norbs, 0.0); 

    while(!f1.eof())
    {
        std::getline(f1, msg);
        boost::algorithm::trim(msg);
        boost::split(tok, msg, boost::is_any_of(" \t"), boost::token_compress_on);
        if (tok.size() == 5)
        {
            double integral = atof(tok[0].c_str());
            int a = atoi(tok[1].c_str()), b = atoi(tok[2].c_str()), c = atoi(tok[3].c_str()), d = atoi(tok[4].c_str());
            if (integral == 0.0 && a == b && b == c && c == d && d == 0)
            {
                break;
            }
            else if (a == b && b == c && c == d && d == 0)
            {
                core_e = integral;
            }
            else if (a != 0 && b != 0 && c == d && d == 0)
            {
                I1.store.at((a - 1) * norbs + (b - 1)) = integral;
                I1.store.at((b - 1) * norbs + (a - 1)) = integral;
            }
            else if (a != 0 && b != 0 && c != 0 && d != 0)
            {
                I2((a - 1), (b - 1), (c - 1), (d - 1)) = integral;
            }

        }
    } // while

    index = 0;
    while(!f2.eof())
    {
        std::getline(f2, msg);
        boost::algorithm::trim(msg);
        boost::split(tok, msg, boost::is_any_of(" \t"), boost::token_compress_on);
        if (tok.size() == norbs + 1)
        {
            for (int i = 0; i < norbs; i++)
            {
                S.store.at(index * norbs + i) = atof(tok[i + 1].c_str());
            }
            index += 1;
            }
    } // while

    f3 >> nuc_e;

    f1.close();
    f2.close();
    f3.close();
} // end readIntegrals
