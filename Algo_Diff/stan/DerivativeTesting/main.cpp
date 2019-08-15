#include <Eigen/Dense>
#include <stan/math.hpp>

using namespace std;
using namespace Eigen;

template <typename T>
T Determinant(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &Matrix)
{
    Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solver(Matrix);
    return solver.determinant();
}

int main()
{
    int n = 10;
    MatrixXd Rd = MatrixXd::Random(n, n);
    double detd = Determinant(Rd);
    cout << detd << endl;

    Matrix<stan::math::var, Dynamic, Dynamic> Rs = Rd;
    stan::math::var dets = Determinant(Rd);
    cout << dets << endl;
}
