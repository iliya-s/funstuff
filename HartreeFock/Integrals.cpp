#include "Integrals.h"
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

void ReadFCIDUMP(std::string FCIDUMP, Integral::OneElectron &I1, Integral::TwoElectron &I2, double &core_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep)
{
    std::ifstream f(FCIDUMP);
    if (!f.is_open())
    {
        std::cout << "Integral file " << FCIDUMP << " does not exist." << std::endl;
        std::exit(1);
    }

    norbs = -1;
    nelec = -1;
    sz = -1;
    I2.ksym == false;

    std::string line;
    std::vector<std::string> token;
    while (!f.eof())
    {
        std::getline(f, line);
        boost::algorithm::trim(line);
        boost::split(token, line, boost::is_any_of(", \t="), boost::token_compress_on);
        if (token.size() == 1 && token[0] == "&END")
        {
            break;
        }
        else
        {
            if (token[0] == "&FCI")
            {
                if (token[1] == "NORB")
                {
                    norbs = std::stoi(token[2]);
                }
                if (token[3] == "NELEC")
                {
                    nelec = std::stoi(token[4]);
                }
                if (token[5] == "MS2")
                {
                    sz = std::stoi(token[6]);
                }
            }
            else if (token[0] == "ISYM")
            {
                continue;
            }
            else if (token[0] == "KSYM")
            {
                I2.ksym = true;
            }
            else if (token[0] == "ORBSYM")
            {
                for (int i = 1; (i < token.size() && token[i] != "\0"); i++)
                {
                    irrep.push_back(std::stoi(token[i]));
                }
            }
            else 
            {
                for (int i = 0; (i < token.size() && token[i] != "\0"); i++)
                {
                    irrep.push_back(std::stoi(token[i]));
                }
            }
        }
    } //while

    if (norbs == -1 || nelec == -1 || sz == -1)
    {
        std::cout << "failed reading norbs, nelec, or sz" << std::endl;
        std::exit(2);
    }

    nalpha = nelec / 2 + sz;
    nbeta = nelec - nalpha;
    irrep.resize(norbs);

    int npair = norbs * (norbs + 1) / 2;
    if (I2.ksym)
    {
        npair = norbs * norbs;
    }

    I2.norbs = norbs;
    I2.store.resize(npair * (npair + 1) / 2);

    I1.norbs = norbs;
    I1.store.resize(norbs * norbs);

    while(!f.eof())
    {
        std::getline(f, line);
        boost::algorithm::trim(line);
        boost::split(token, line, boost::is_any_of(" \t"), boost::token_compress_on);
        if (token.size() == 5)
        {
            double integral = std::stod(token[0]);
            int a = std::stoi(token[1]);
            int b = std::stoi(token[2]);
            int c = std::stoi(token[3]);
            int d = std::stoi(token[4]);
            if (a == 0 && b == 0 && c == 0 && d == 0)
            {
                core_e = integral;
            }
            else if (b == 0 && c == 0 && d == 0)
            {
                continue;   //orbital energy
            }
            else if (a != 0 && b != 0 && c == 0 && d == 0)
            {
                I1((a - 1), (b - 1)) = integral;
                I1((b - 1), (a - 1)) = integral;
            }
            else if (a != 0 && b != 0 && c != 0 && d != 0)
            {
                I2((a - 1), (b - 1), (c - 1), (d - 1)) = integral;
            }
        }
    } // while

    I2.Direct.resize(norbs, norbs);
    I2.Exchange.resize(norbs, norbs);
    for (int i = 0; i < norbs; i++)
    {
        for (int j = 0; j < norbs; j++)
        {
            I2.Direct(i, j) = I2(i, i, j, j);
            I2.Exchange(i, j) = I2(i, j, j, i);
        }
    }
    f.close();
}

void ReadSquareMatrixIntoOneInt(std::string MatrixFile, int norbs, Integral::OneElectron &I)
{
    std::ifstream f(MatrixFile);
    if (!f.is_open())
    {
        std::cout << "Matrix file " << MatrixFile << " does not exist." << std::endl;
        exit(1);
    }

    I.norbs = norbs;
    I.store.resize(norbs * norbs);
    for (int i = 0; i < norbs; i++)
    {
        for (int j = 0; j < norbs; j++)
        {
            f >> I(i, j);
        }
    }

    f.close();
}

void ReadMatrix(std::string MatrixFile, int rows, int cols, Eigen::MatrixXd &M)
{
    std::ifstream f(MatrixFile);
    if (!f.is_open())
    {
        std::cout << "Matrix file " << MatrixFile << " does not exist." << std::endl;
        exit(1);
    }

    M.resize(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            f >> M(i, j);
        }
    }

    f.close();
}

void WriteMatrix(std::string MatrixFile, const Eigen::MatrixXd &M)
{
    std::ofstream f(MatrixFile);
    if (!f.is_open())
    {
        std::cout << "Matrix file " << MatrixFile << " could not be opened" << std::endl;
    }
    for (int i = 0; i < M.rows(); i++)
    {
        for (int j = 0; j < M.cols(); j++)
        {
            f << M(i, j) << "\t";
        }
        f << std::endl;
    }
}

void ReadAtomicOrbitalIntegrals(std::string AOFCIDUMP, std::string METRIC, Integral::OneElectron &I1, Integral::TwoElectron &I2, Integral::OneElectron &S, double &core_e, int &norbs, int &nelec, int &nalpha, int &nbeta, int &sz, std::vector<int> &irrep)
{
    ReadFCIDUMP(AOFCIDUMP, I1, I2, core_e, norbs, nelec, nalpha, nbeta, sz, irrep);
    ReadSquareMatrixIntoOneInt(METRIC, norbs, S);
}
