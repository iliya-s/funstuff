#include <vector>
#include <iostream>
#include "MergeSort.h"
#include "BubbleSort.h"
#include "InsertionSort.h"

template<typename T>
class NumericComparison
{
    public:
    bool operator()(T LHS, T RHS)
    {
        if (LHS < RHS)
            return true;
        else
            return false;
    }
};

int main(int agrc, char **argv)
{
    /*
    std::vector<int> a = { 10};
    std::vector<int> b = { 2, 8, 11, 12, 13};
    std::vector<int> c;
    MergeSortedVectors(a, b, c);
    for (int i = 0; i < c.size(); i++)
    {
        std::cout << c[i] << " ";
    }
    std::cout << std::endl;
    */
    
    std::vector<int> test = { 10, 4, 2, 8, 6, 11, 155, 123, 67, 89, 34, 0 , -2, 3, 8, 100, 43, 66, 20 };
    std::cout << "Array" << std::endl;
    for (int i = 0; i < test.size(); i++)
    {
        std::cout << test[i] << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "MergeSort" << std::endl;
    MergeSort(test);
    for (int i = 0; i < test.size(); i++)
    {
        std::cout << test[i] << " ";
    }
    std::cout << std::endl << std::endl;
    

    std::cout << "BubbleSort" << std::endl;
    std::vector<int> test1 = { 10, 4, 2, 8, 6, 11, 155, 123, 67, 89, 34, 0 , -2, 3, 8, 100, 43, 66, 20 };
    BubbleSort(test1);
    for (int i = 0; i < test1.size(); i++)
    {
        std::cout << test1[i] << " ";
    }
    std::cout << std::endl << std::endl;
    

    std::cout << "InsertionSort" << std::endl;
    std::vector<int> test2 = { 10, 4, 2, 8, 6, 11, 155, 123, 67, 89, 34, 0 , -2, 3, 8, 100, 43, 66, 20 };
    InsertionSort(test2);
    for (int i = 0; i < test2.size(); i++)
    {
        std::cout << test2[i] << " ";
    }
    std::cout << std::endl << std::endl;
}
