#ifndef BUBBLE_SORT_HEADER
#define BUBBLE_SORT_HEADER
#include <vector>

template<typename T>
void BubbleSort(std::vector<T> &V)
{
    while(true)
    {
        int swaps = 0;
        for (int i = 0; i < (V.size() - 1); i++)
        {
            if (V[i + 1] < V[i])
            {
                T buffer = V[i];
                V[i] = V[i + 1];
                V[i + 1] = buffer;
                swaps++;
            }
            else
                continue;
        }
        if (swaps == 0)
            break;
    }
}
#endif
