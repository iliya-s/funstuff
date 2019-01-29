#ifndef MERGE_SORT_HEADER_H
#define MERGE_SORT_HEADER_H
#include <vector>

template<typename T>
void MergeSortedVectors(const std::vector<T> &LHS, const std::vector<T> &RHS, std::vector<T> &OUT)
{
    OUT.clear();
    int l = 0, r = 0;
    while (l < LHS.size() && r < RHS.size())
    {
        if (LHS[l] < RHS[r])
        {
            OUT.push_back(LHS[l]);
            l++;
        }
        else
        {
            OUT.push_back(RHS[r]);
            r++;
        }
    }
    while (l < LHS.size())
    {
        OUT.push_back(LHS[l]);
        l++;
    }
    while (r < RHS.size())
    {
        OUT.push_back(RHS[r]);
        r++;
    }
}

template<typename T>
void MergeSort(std::vector<T> &V)
{
    if (V.size() == 1)
    {
        return;
    }
    else 
    {
        int mid = V.size() / 2;
        std::vector<T> LHS(V.begin(), V.begin() + mid);
        std::vector<T> RHS(V.begin() + mid, V.end());
        MergeSort(LHS);
        MergeSort(RHS);
        MergeSortedVectors(LHS, RHS, V);
    }
}
#endif
