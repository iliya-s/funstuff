#ifndef BUBBLESORT_HEADER_H
#define BUBBLESORT_HEADER_H
#include <vector>

template<typename T>
void BubbleSort(std::vector<T> &V)
{
    while(true)
    {
        int swaps = 0;
        for (int i = 0; i < V.size() - 1; i++)
        {
            if (V.at(i + 1) < V.at(i))
            {
                std::swap(V.at(i), V.at(i + 1));
                swaps++;
            }
        }
        if (swaps == 0)
        {
            break;
        }
    }
}
#endif
