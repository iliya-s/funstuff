#ifndef INSERTION_SORT_HEADER_H
#define INSERTION_SORT_HEADER_H
#include <vector>

template<typename T>
void InsertionSort(std::vector<T> &V)
{
    for (int i = 1; i < V.size(); i++)
    {
        T val = V[i];
        int j = i - 1;
        while (j >= 0 && V[j] > val)
        {
            V[j + 1] = V[j];
            j--;
        }
        V[j + 1] = val;
    }
}
#endif
