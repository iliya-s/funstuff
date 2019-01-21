#ifndef HAMILTONIAN_OPERATORS_HEADER
#define HAMILTONIAN_OPERATORS_HEADER
#include "Determinant.h"

//calculates the number of different occupied orbitals between two determinants
int NumDiffOrbs(const Determinant &LHS, const Determinant &RHS)
{
    int count = 0;
    for (int i = 0; i < RHS.len; i++)
    {
        count += CountSetBits(LHS.String[0][i] ^ RHS.String[0][i]);
        count += CountSetBits(LHS.String[1][i] ^ RHS.String[1][i]);
    }
    return count;
}

//finds spin orbital indices of differing occupied orbitals for two determinants with 2 different orbitals occupied
//  i corresponds to the orbital occupied in LHS, a corresponds to the orbital occupied in RHS
void DiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &a)
{
    int indices[2];
    for (int l = 0; l < LHS.len; l++)
    {
        long alpha = LHS.String[0][l] ^ RHS.String[0][l];
        long beta = LHS.String[1][l] ^ RHS.String[1][l];
        int pos = 0;
        while (alpha != 0 && beta != 0)
        {
            indices[pos] = 2.0 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
            indices[pos] = 2.0 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
        while (alpha != 0)
        {
            indices[pos] = 2.0 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
        }
        while (beta != 0)
        {
            indices[pos] = 2.0 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
    }
    /*
    int indices[2];
    int spin[2];
    long test[2];
    for (int l = 0; l < RHS.len; l++)
    {
        long alpha = LHS.String[0][l] ^ RHS.String[0][l];
        long beta = LHS.String[1][l] ^ RHS.String[1][l];
        if (alpha == 0 && beta == 0) //no difference between bitstring
        {
            continue;
        }
        else if (alpha != 0 && beta == 0)   //alpha excitation
        {
            test[0] = alpha;
            test[1] = alpha ^ (alpha & -alpha);
            spin[0] = 0;
            spin[1] = 0;
        }
        else if (alpha == 0 && beta != 0)   //beta excitation
        {
            test[0] == beta;
            test[1] = beta ^ (beta & -beta);
            spin[0] = 1;
            spin[1] = 1;
        }
        else if (alpha != 0 && beta != 0)   //mixed spin excitation
        {
            test[0] = alpha & -alpha;
            test[1] = beta & -beta;
            spin[0] = 0;
            spin[1] = 1;
        }
        for (int m = 0; m < 2; m++)
        {
            if (test[m] != 0)
                indices[m] = 2.0 * (__builtin_ctzl(test[m]) + l * 64) + spin[m];
        }
    }
    */
    if (LHS(indices[0])) //LHS has ith orbital occupied
    {
        i = indices[0];
        a = indices[1];
    }
    else
    {
        i = indices[1];
        a = indices[0];
    }   
}

//finds spin orbital indices of differing occupied orbitals for two determinants with 4 different orbitals occupied
//  i, j corresponds to the orbitals occupied in LHS, a, b corresponds to the orbital occupied in RHS
void DiffOrbIndices(const Determinant &LHS, const Determinant &RHS, int &i, int &j, int &a, int &b)
{
    int indices[4];
    for (int l = 0; l < LHS.len; l++)
    {
        long alpha = LHS.String[0][l] ^ RHS.String[0][l];
        long beta = LHS.String[1][l] ^ RHS.String[1][l];
        int pos = 0;
        while (alpha != 0 && beta != 0)
        {
            indices[pos] = 2.0 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
            indices[pos] = 2.0 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
        while (alpha != 0)
        {
            indices[pos] = 2.0 * (__builtin_ctzl(alpha) + l * 64) + 0;
            long ra = alpha & -alpha;
            alpha ^= ra;
            pos++;
        }
        while (beta != 0)
        {
            indices[pos] = 2.0 * (__builtin_ctzl(beta) + l * 64) + 1;
            long rb = beta & -beta;
            beta ^= rb;
            pos++;
        }
    }
    if (LHS(indices[0]) && LHS(indices[1]))
    {
        i = std::min(indices[0], indices[1]);
        j = std::max(indices[0], indices[1]);
        a = std::min(indices[2], indices[3]);
        b = std::max(indices[2], indices[3]);
    }
    else if (LHS(indices[0]) && LHS(indices[2]))
    {
        i = std::min(indices[0], indices[2]);
        j = std::max(indices[0], indices[2]);
        a = std::min(indices[1], indices[3]);
        b = std::max(indices[1], indices[3]);
    }
    else if (LHS(indices[0]) && LHS(indices[3]))
    {
        i = std::min(indices[0], indices[3]);
        j = std::max(indices[0], indices[3]);
        a = std::min(indices[1], indices[2]);
        b = std::max(indices[1], indices[2]);
    }
    else if (LHS(indices[1]) && LHS(indices[2]))
    {
        i = std::min(indices[1], indices[2]);
        j = std::max(indices[1], indices[2]);
        a = std::min(indices[0], indices[3]);
        b = std::max(indices[0], indices[3]);
    }
    else if (LHS(indices[1]) && LHS(indices[3]))
    {
        i = std::min(indices[1], indices[3]);
        j = std::max(indices[1], indices[3]);
        a = std::min(indices[2], indices[0]);
        b = std::max(indices[2], indices[0]);
    }
    else if (LHS(indices[2]) && LHS(indices[3]))
    {
        i = std::min(indices[2], indices[3]);
        j = std::max(indices[2], indices[3]);
        a = std::min(indices[1], indices[0]);
        b = std::max(indices[1], indices[0]);
    }
}


#endif
