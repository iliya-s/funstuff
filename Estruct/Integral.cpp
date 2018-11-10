#include "Integral.h"
#include <fstream>
#include <math.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string.hpp>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

using namespace std;
using namespace boost;

int read_ao_integrals(string ao_fcidump, string metric_file, string nuc_file, OneInt &AI1, TwoInt &AI2, OneInt &S, double &core_e, double &nuc_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, vector<int> &irrep)
{
    ifstream f1(ao_fcidump);
    if (!f1.is_open())
    {
        cout << "Integral file " << ao_fcidump << " does not exist." << endl;
        return 1;
    }
    ifstream f2(metric_file);
    if (!f2.is_open())
    {
        cout << "Integral file " << metric_file << " does not exist." << endl;
        return 2;
    }
    ifstream f3(nuc_file);
    if (!f3.is_open())
    {
        cout << "Nuclear energy file " << nuc_file << " does not exist." << endl;
        return 3;
    }

    bool start_scaling = false;
    norbs = -1;
    nelec = -1;
    sz = -1;

    int index = 0;
    vector<string> token;
    string msg;
    while (!f1.eof())
    {
        getline(f1, msg);
        algorithm::trim(msg);
        split(token, msg, is_any_of(", \t="), token_compress_on);
        if (start_scaling == false && token.size() == 1 && (iequals(token[0], "&END") || iequals(token[0], "/")))
        {
            start_scaling = true;
            index += 1;
            break;
        }
        else if (start_scaling == false)
        {
            if (iequals(token[0].substr(0,4), "&FCI"))
            {
                if (iequals(token[1].substr(0,4), "NORB"))
                {
                    norbs = atoi(token[2].c_str());
                }
                if (iequals(token[3].substr(0,5), "NELEC"))
                {
                    nelec = atoi(token[4].c_str());
                }
                if (iequals(token[5].substr(0,5), "MS2"))
                {
                    sz = atoi(token[6].c_str());
                }
            } 
            else if (iequals(token[0].substr(0,4), "ISYM"))
            {
                continue;
            }
            else if (iequals(token[0].substr(0,4), "KSYM"))
            {
                continue;
            }
            else if (iequals(token[0].substr(0,6), "ORBSYM"))
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
        cout << "failed reading norbs, nelec, or sz" << endl;
        return 4;
    }
    nalpha = nelec / 2 + sz;
    nbeta = nelec - nalpha;
    irrep.resize(norbs);
    core_e = 0.0;

    AI2.store.clear();
    AI2.norbs = norbs;
    AI2.store.resize((norbs * norbs * norbs * norbs), 0.0);

    AI1.norbs = norbs;
    AI1.store.setZero(norbs,norbs);

    S.store.setZero(norbs,norbs);
    S.norbs = norbs;

    while(!f1.eof())
    {
        getline(f1, msg);
        algorithm::trim(msg);
        split(token, msg, is_any_of(" \t"), token_compress_on);
        if (token.size() == 5)
        {
            double integral = atof(token[0].c_str());
            int a = atoi(token[1].c_str()), b = atoi(token[2].c_str()), c = atoi(token[3].c_str()), d = atoi(token[4].c_str());
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
                AI1((a - 1), (b - 1)) = integral;
                AI1((b - 1), (a - 1)) = integral;
            }
            else if (a != 0 && b != 0 && c != 0 && d != 0)
            {
                AI2((a - 1), (b - 1), (c - 1), (d - 1)) = integral;
            }

        }
    } // while

    index = 0;
    while(!f2.eof())
    {
        getline(f2, msg);
        algorithm::trim(msg);
        split(token, msg, is_any_of(" \t"), token_compress_on);
        if (token.size() > norbs)
        {
            for (int i = 0; i < norbs; i++)
            {
                S(index, i) = atof(token[i + 1].c_str());
            }
            index += 1;
        }
    } // while

    f3 >> nuc_e;

    f1.close();
    f2.close();
    f3.close();
    return 0;
} // end readIntegrals
