#include <iostream>
#include <random>
#include <time.h>

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
    double count = 0.0;
    for (int i = 0; i < niter; ++i)
    {
        double x = dist(rng);
        double y = dist(rng);
        double test = x * x + y * y;
        if (test <= 1.0) count += 1.0;
    }
    double pi = 4.0 * count / (double) niter;

#ifndef SERIAL
    double finalpi;
    boost::mpi::reduce(world, pi, finalpi, std::plus<double> (), 0);
    finalpi /= world.size();
    if (world.rank() == 0)
    {
        std::cout << "estimate for pi: " << finalpi << std::endl;
    }
#else
    std::cout << "estimate for pi: " << pi << std::endl;
#endif
    return 0;
}
