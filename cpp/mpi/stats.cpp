#include "stats.h"

template <class T>
double mean(std::vector<T> v)
{
    double sum = 0.0;
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        sum += (double)*it;
    }
    double avg = sum / (double) v.size();
    return avg;
}

template <class T>
double variance(std::vector<T> v, double avg)
{
    double ss = 0.0;
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        ss += (*it - avg) * (*it - avg);
    }
    double var = ss / (double) v.size();
    return var;
}

template
double mean<double> (std::vector<double> v);
template
double variance<double> (std::vector<double> v, double avg);

template
double mean<float> (std::vector<float> v);
template
double variance<float> (std::vector<float> v, double avg);

template
double mean<int> (std::vector<int> v);
template
double variance<int> (std::vector<int> v, double avg);
