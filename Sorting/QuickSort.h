#ifndef QUICKSORT_HEADER_H
#define QUICKSORT_HEADER_H
#include <vector>
#include <random>

template <typename T>
int Partition(std::vector<T> &V, int l, int r)
{
    /*
    T pivot = V.at(r);
    int i = l;
    for (int j = l; j < r; j++)
    {
        if (V.at(j) < pivot)
        {
            std::swap(V.at(i), V.at(j));
            i++;
        }
    }
    std::swap(V.at(i), V.at(r));
    return i;
    */
    T pivot = V.at(r);
    int i = 0;
    int j = 0;
    while (j < r)
    {
        if (V.at(j) <= pivot)
        {
            std::swap(V.at(j), V.at(i));
            i++;
            j++;
        }
        else if (V.at(j) > pivot)
        {
            j++;
        }
    }
    std::swap(V.at(r), V.at(i));
    return i;
}

template <typename T>
int RandomPartition(std::vector<T> &V, int l, int r)
{
    std::mt19937 rng(1);
    std::uniform_int_distribution<int> dist(l, r);
    std::swap(V.at(dist(rng)), V.at(r));
    return Partition(V, l, r);
}

template <typename T>
void QuickSortHelper(std::vector<T> &V, int l, int r)
{
    if (l < r)
    {
        int p = Partition(V, l, r);
        QuickSortHelper(V, l, p - 1);
        QuickSortHelper(V, p + 1, r);
    }
}

template <typename T>
void RandomQuickSortHelper(std::vector<T> &V, int l, int r)
{
    if (l < r)
    {
        int p = RandomPartition(V, l, r);
        QuickSortHelper(V, l, p - 1);
        QuickSortHelper(V, p + 1, r);
    }
}

template <typename T>
void QuickSort(std::vector<T> &V)
{
    QuickSortHelper(V, 0, V.size() - 1);
}

template <typename T>
void RandomQuickSort(std::vector<T> &V)
{
    QuickSortHelper(V, 0, V.size() - 1);
}
#endif
