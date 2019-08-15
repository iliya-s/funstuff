#ifndef INSERTIONSORT_HEADER_H
#define INSERTIONSORT_HEADER_H
#include <vector>

template<typename T>
void InsertionSort(std::vector<T> &V)
{
    for (int i = 1; i < V.size(); i++)
    {
        T val = V.at(i);
        int j = i - 1;
        while (j >= 0 && V.at(j) > val)
        {
            V.at(j + 1) = V.at(j);
            j--;
        }
        V.at(j + 1) = val;
    }
}
#endif
