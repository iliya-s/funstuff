#include <iostream>
#include <random>
#include <time.h>
#include <vector>
#include "stats.h"

#ifndef SERIAL
#include <boost/mpi.hpp>
#endif

int main(int argc, char *argv[])
{
    std::mt19937 rng;
#ifndef SERIAL
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    rng.seed(time(nullptr) + world.rank());
#else
    rng.seed(time(nullptr));
#endif
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    int niter = atoi(argv[1]);
    std::vector<double> v(niter);
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        double Sn = 0.0;
        double n = 0.0;
        while (true)
        {
            Sn += dist(rng);
            n += 1.0;
            if (Sn > 1.0) break;
        }
        *it = n;
    }
    double avg = mean<double> (v);
    double var = variance<double> (v,avg);

#ifndef SERIAL
    double finalavg;
    boost::mpi::reduce(world, avg, finalavg, std::plus<double> (), 0);
    finalavg /= world.size();
    if (world.rank() == 0)
    {
        std::cout << "estimate for e: " << finalavg << std::endl;
    }
#else
    std::cout << "estimate for e: " << avg << std::endl;
#endif
    return 0;
}
