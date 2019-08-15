#include <Eigen/Dense>
#include <stan/math.hpp>
#include <complex>
#include "Complex.h"

using namespace std;
using namespace Eigen;

int main()
{
    cout << "__std::complex__" << endl;
    std::complex<double> a0;
    cout << a0 << endl;
    a0.imag(1.0);
    a0.real(1.0);
    std::cout << a0 << std::endl;
    std::complex<double> a1(2.0);
    std::complex<double> a2(2.0,3.0);
    cout << a1 << endl;
    cout << std::conj(a2) << endl;
    cout << a2 / (a0 * a1) + a0 - a1 << endl;
    cout << endl;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A = MatrixXcd::Random(10, 10), A1(10, 10);
    A1 = A.adjoint();
    //cout << (A * A1).real() << endl;

    Eigen::FullPivLU<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> lua(A);
    A1 = lua.inverse();
    std::complex<double> deta = lua.determinant();
    cout << deta << endl << endl;
    cout << A1.real() << endl;
    cout << endl << A1.imag() << endl;

    Eigen::HessenbergDecomposition<Eigen::MatrixXcd> hda(A);
    cout << endl << hda.matrixH().real() << endl;

    cout << "__Complex__" << endl;
    cout << "double" << endl;
    Complex<double> b0;
    std::cout << b0 << std::endl;
    b0.imag(1.0);
    b0.real(1.0);
    std::cout << b0 << std::endl;
    Complex<double> b1(2.0);
    Complex<double> b2(2.0,3.0);
    cout << b1 << endl;
    cout << b2.conj() << endl;
    cout << b2 / (b0 * b1) + b0 - b1 << endl;
    cout << endl;

    Eigen::Matrix<Complex<double>, Eigen::Dynamic, Eigen::Dynamic> B(10, 10), B1(10, 10);
    for (int i = 0; i < B.rows(); i++)
    {
        for (int j = 0; j < B.cols(); j++)
        {
            B(i, j) = A(i, j);
        }
    }
    B1 = B.adjoint();
    //cout << (B * B1).real() << endl;

    Eigen::FullPivLU<Eigen::Matrix<Complex<double>, Eigen::Dynamic, Eigen::Dynamic>> lub(B);
    B1 = lub.inverse();
    Complex<double> detb = lub.determinant();
    cout << detb << endl << endl;
    cout << B1.real() << endl;
    cout << endl << B1.imag() << endl;

    Eigen::HessenbergDecomposition<Eigen::Matrix<Complex<double>, Eigen::Dynamic, Eigen::Dynamic>> hdb(B);
    cout << endl << hdb.matrixH().real() << endl;

    cout << "stan::math::var" << endl;
    Complex<stan::math::var> c0;
    std::cout << c0 << std::endl;
    c0.imag(1.0);
    c0.real(1.0);
    std::cout << c0 << std::endl;
    Complex<stan::math::var> c1(2.0);
    Complex<stan::math::var> c2(2.0,3.0);
    cout << c1 << endl;
    cout << c2.conj() << endl;
    cout << c2 / (c0 * c1) + c0 - c1 << endl;
    cout << endl;

    Eigen::Matrix<Complex<stan::math::var>, Eigen::Dynamic, Eigen::Dynamic> C(10, 10), C1(10, 10);
    for (int i = 0; i < C.rows(); i++)
    {
        for (int j = 0; j < C.cols(); j++)
        {
            C(i, j) = A(i, j);
        }
    }
    C1 = C.adjoint();
    //cout << (C * C1).real() << endl;

    Eigen::FullPivLU<Eigen::Matrix<Complex<stan::math::var>, Eigen::Dynamic, Eigen::Dynamic>> luc(C);
    C1 = luc.inverse();
    Complex<stan::math::var> detc = luc.determinant();
    cout << detc << endl << endl;
    cout << C1.real() << endl;
    cout << endl << C1.imag() << endl;

    Eigen::HessenbergDecomposition<Eigen::Matrix<Complex<stan::math::var>, Eigen::Dynamic, Eigen::Dynamic>> hdc(C);
    cout << endl << hdc.matrixH().real() << endl;
}
