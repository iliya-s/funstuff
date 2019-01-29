#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;

void conj_grad(MatrixXd &A, VectorXd &b, VectorXd &xguess, int n, VectorXd &x)
{
    VectorXd x0 = xguess;
    VectorXd r0 = b - A * x0;
    VectorXd p0 = r0;
    double r0_square = r0.adjoint() * r0;
    for (int i = 0; i < n; i++)
    {
        VectorXd Ap = A * p0;
        double pAp = p0.adjoint() * Ap;
        double alpha = r0_square / pAp;
        VectorXd x1 = x0 + alpha * p0;
        VectorXd r1 = r0 - alpha * Ap;
        double r1_square = r1.adjoint() * r1;
        double beta = r1_square / r0_square;
        VectorXd p1 = r1 + beta * p0;

        x0 = x1;
        r0 = r1;
        p0 = p1;
        r0_square = r1_square;
        x = x1;
    }
}

void pinv(MatrixXd &A, MatrixXd &Ainv)
{
    SelfAdjointEigenSolver<MatrixXd> es;
    es.compute(A);
    std::vector<int> cols;
    for (int m = 0; m < A.cols(); m++)
    {
      if (fabs(es.eigenvalues()(m)) > 1.e-9)
      {
        cols.push_back(m);
      }
    }
    MatrixXd U = MatrixXd::Zero(A.rows(), cols.size());
    MatrixXd eig_inv = MatrixXd::Zero(cols.size(),cols.size());
    for (int m = 0; m < cols.size(); m++)
    {
      int index = cols[m];
      U.col(m) = es.eigenvectors().col(index);
      double eigval = es.eigenvalues()(index);
      eig_inv(m,m) = 1.0 / eigval;
    }
    Ainv = U * eig_inv * U.adjoint();
}

int main()
{
    int D = 5;
    MatrixXd A = MatrixXd::Random(D,D);
    MatrixXd B = A.adjoint();
    A = A + B;
    VectorXd b = VectorXd::Random(D);

    MatrixXd Ainv = MatrixXd::Zero(D,D);
    pinv(A, Ainv);
    std::cout << Ainv * A << std::endl;
    std::cout << Ainv * b << std::endl;
    
    VectorXd xguess = VectorXd::Random(D);
    VectorXd x = VectorXd::Zero(D);
    conj_grad(A, b, xguess, 10, x);
    std::cout << x << std::endl;

}
