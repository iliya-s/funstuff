#include <Eigen/Dense>
#include <stan/math.hpp>
#include <complex>

int main()
{
    std::complex<double> x(1,1);
    auto a = 2.0 * x;
    std::cout << a << std::endl;
}
