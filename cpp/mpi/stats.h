#ifndef STATS
#define STATS

#include <vector>

template <class T>
double mean(std::vector<T> v);

template <class T>
double variance(std::vector<T> v, double avg);

#endif
