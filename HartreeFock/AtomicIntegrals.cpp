#include "AtomicIntegrals.h"
#include <fstream>
#include <math.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string.hpp>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

int ReadAOIntegrals(std::string AOfcidump, std::string metric_file, std::string nuc_file, AOneInt &AI1, ATwoInt &AI2, AOneInt &S, double &core_e, double &nuc_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep)
{
    std::ifstream f1(AOfcidump);
    if (!f1.is_open())
    {
        std::cout << "Integral file " << AOfcidump << " does not exist." << std::endl;
        return 1;
    }
    std::ifstream f2(metric_file);
    if (!f2.is_open())
    {
        std::cout << "Integral file " << metric_file << " does not exist." << std::endl;
        return 2;
    }
    std::ifstream f3(nuc_file);
    if (!f3.is_open())
    {
        std::cout << "Nuclear energy file " << nuc_file << " does not exist." << std::endl;
        return 3;
    }

    bool start_scaling = false;
    norbs = -1;
    nelec = -1;
    sz = -1;

    int index = 0;
    std::vector<std::string> token;
    std::string msg;
    while (!f1.eof())
    {
        std::getline(f1, msg);
        boost::algorithm::trim(msg);
        boost::split(token, msg, boost::is_any_of(", \t="), boost::token_compress_on);
        if (start_scaling == false && token.size() == 1 && (boost::iequals(token[0], "&END") || boost::iequals(token[0], "/")))
        {
            start_scaling = true;
            index += 1;
            break;
        }
        else if (start_scaling == false)
        {
            if (boost::iequals(token[0].substr(0,4), "&FCI"))
            {
                if (boost::iequals(token[1].substr(0,4), "NORB"))
                    norbs = atoi(token[2].c_str());
                if (boost::iequals(token[3].substr(0,5), "NELEC"))
                    nelec = atoi(token[4].c_str());
                if (boost::iequals(token[5].substr(0,5), "MS2"))
                    sz = atoi(token[6].c_str());
            }
            else if (boost::iequals(token[0].substr(0,4), "ISYM"))
                continue;
            else if (boost::iequals(token[0].substr(0,4), "KSYM"))
                continue;
            else if (boost::iequals(token[0].substr(0,6), "ORBSYM"))
            {
                for (int i = 1; i < token.size(); i++)
                {
                    irrep.push_back(atoi(token[i].c_str()));
                }
            }
            else 
            {
                for (int i = 0; i < token.size(); i++)
                {
                    irrep.push_back(atoi(token[i].c_str()));
                }
            }
            index += 1;
        }
    } //while

    if (norbs == -1 || nelec == -1 || sz == -1)
    {
        std::cout << "failed reading norbs, nelec, or sz" << std::endl;
        return 4;
    }
    nalpha = nelec / 2 + sz;
    nbeta = nelec - nalpha;
    irrep.resize(norbs);
    core_e = 0.0;

    AI2.store.clear();
    AI2.norbs = norbs;
    AI2.store.assign(norbs * norbs * norbs * norbs, 0.0);

    AI1.norbs = norbs;
    AI1.store.assign(norbs * norbs, 0.0);

    S.store.assign(norbs * norbs, 0.0);
    S.norbs = norbs;

    while(!f1.eof())
    {
        std::getline(f1, msg);
        boost::algorithm::trim(msg);
        boost::split(token, msg, boost::is_any_of(" \t"), boost::token_compress_on);
        if (token.size() == 5)
        {
            double integral = atof(token[0].c_str());
            int a = atoi(token[1].c_str()), b = atoi(token[2].c_str());
            int c = atoi(token[3].c_str()), d = atoi(token[4].c_str());
            if (integral == 0.0 && a == 0 && b == 0 && c == 0 && d == 0)
                break;
            else if (a == 0 && b == 0 && c == 0 && d == 0)
                core_e = integral;
            else if (a != 0 && b != 0 && c == 0 && d == 0)
            {
                AI1((a - 1), (b - 1)) = integral;
                AI1((b - 1), (a - 1)) = integral;
            }
            else if (a != 0 && b != 0 && c != 0 && d != 0)
                AI2((a - 1), (b - 1), (c - 1), (d - 1)) = integral;

        }
    } // while

    index = 0;
    while(!f2.eof())
    {
        std::getline(f2, msg);
        boost::algorithm::trim(msg);
        boost::split(token, msg, boost::is_any_of(" \t"), boost::token_compress_on);
        if (token.size() > norbs)
        {
            for (int i = 0; i < norbs; i++)
                S(index, i) = atof(token[i + 1].c_str());
            index += 1;
        }
    } // while

    f3 >> nuc_e;

    f1.close();
    f2.close();
    f3.close();
    return 0;
} // end ReadAOIntegrals
