#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

using namespace Eigen;
using namespace std;

static std::mt19937 rg(1);

//sorts the eigenvalues in D and their respective eigenvectors in V in increasing order
void SortEig(Eigen::VectorXcd &D, Eigen::MatrixXcd &V)
{
  Eigen::VectorXcd Dcopy(D);
  Eigen::MatrixXcd Vcopy(V);
  std::vector<std::complex<double>> a(Dcopy.data(), Dcopy.data() + Dcopy.size());
  // initialize original index locations
  std::vector<int> idx(a.size());
  std::iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in a
  //std::sort(idx.begin(), idx.end(), [&a](int i1, int i2) -> bool { return std::abs(a[i1]) < std::abs(a[i2]); });
  std::sort(idx.begin(), idx.end(), [&a](int i1, int i2) -> bool { return a[i1].real() < a[i2].real(); });
  for (int i = 0; i < idx.size(); i++)
  {
    D(i) = Dcopy(idx[i]);
    V.col(i) = Vcopy.col(idx[i]);
  }
}

Eigen::VectorXd RandomVector(int n)
{
    Eigen::VectorXd o(n);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < o.size(); i++) { o(i) = dist(rg); }
    return o;
}

//DOI. 10.1137/090771806
int RandomizedRangeFinder(const Eigen::MatrixXd &AO, Eigen::MatrixXd &Q, double epsilon, int r)
{
    int dim = AO.rows();

    Eigen::MatrixXd Y = AO.block(0, 0, dim, r);
    Eigen::VectorXd y(r);
    for (int i = 0; i < r; i++) { y(i) = Y.col(i).norm(); }

    double tol = epsilon / (10 * std::sqrt(2 / M_PI));
    int j = 0;
    Q = MatrixXd::Zero(dim, 0);
    while (y.maxCoeff() > tol)
    {
        j++;

        Y.col(j - 1) = Y.col(j - 1) - Q * (Q.transpose() * Y.col(j - 1));
        Eigen::VectorXd q = Y.col(j - 1) / (Y.col(j - 1).norm() + 1.0e-8);

        Q.conservativeResize(dim, j);
        Q.col(j - 1) = q;

        Y.conservativeResize(dim, j + r);
        Y.col(j + r - 1) = AO.col(j + r - 1) - Q * (Q.transpose() * AO.col(j + r - 1));

        for (int i = j + 1; i <= j + r - 1; i++)
        {
            Y.col(i - 1) = Y.col(i - 1) - q * (q.transpose() * Y.col(i - 1));
        }

        for (int i = j + 1; i <= j + r; i++)
        {
            y(i - j - 1) = Y.col(i - 1).norm();
        }

        //cout << j << " | " << y.transpose() << endl;
        if (j >= AO.cols() - r) break;
    }
    return j;
}


int main(int argc, char *argv[])
{
    MatrixXd H, S;
    ifstream Hin("H.txt");
    ifstream Sin("S.txt");
    int Hrows, Hcols, Srows, Scols;
    Hin >> Hrows;
    Hin >> Hcols;
    Sin >> Srows;
    Sin >> Scols;
    int dim = Hrows;
    cout << Hrows << " " << Hcols << endl;
    H.resize(Hrows, Hcols);
    S.resize(Srows, Scols);
    for (int i = 0; i < H.rows(); i++)
    {
        for (int j = 0; j < H.cols(); j++)
        {
            Hin >> H(i, j);
            Sin >> S(i, j);
        }
    }
    Hin.close();
    Sin.close();

    cout << "read H, S matrices" << endl << endl;
    //cout << H << endl << endl;
    //cout << S << endl << endl;
    
    auto beginEig = std::chrono::steady_clock::now();

    MatrixXd Sprime = S + 0.00001 * MatrixXd::Identity(dim, dim);
    //MatrixXd Sprime = S;
    MatrixXd Hprime = H + 0.1 * MatrixXd::Identity(dim, dim);
    Eigen::GeneralizedEigenSolver<MatrixXd> es(H, Sprime);

    cout << "Lowest eigenpair of H, S" << endl << endl;
    VectorXcd D = es.eigenvalues();
    MatrixXcd U = es.eigenvectors();
    SortEig(D, U);
    std::complex<double> theta = D(0);
    VectorXcd lowvec = U.col(0);
    cout << theta << endl << endl;
    //cout << lowvec.transpose() << endl << endl;

    cout << "residual" << endl;
    cout << ((H - theta * S) * lowvec).norm() << endl << endl;
    cout << "all eigenvalues" << endl;
    for (int i = 0; i < 10; i++) { cout << D(i) << " "; }
    cout << endl << endl;
    //cout << D.transpose() << endl << endl;

    auto endEig = std::chrono::steady_clock::now();
    std::cout << "Time for EigenSolver: " << std::chrono::duration_cast<std::chrono::seconds>(endEig - beginEig).count() << "[s]" << std::endl;
    
    cout << "Using Random algorithm" << endl << endl;
    //cout << argv[1] << endl;
    auto beginRand = std::chrono::steady_clock::now();

    Eigen::MatrixXd O(dim, 500);
    for (int i = 0; i < O.cols(); i++) { O.col(i) = RandomVector(dim); }
    //double shift = std::stod(argv[1]);
    //Eigen::MatrixXd SO = S * O + 0.00001 * O;
    Eigen::MatrixXd SO = S * O;
    Eigen::MatrixXd Q;
    //double tol = std::stod(argv[1]);
    cout << RandomizedRangeFinder(SO, Q, 0.01, 10) << endl;
    cout << "Q" << endl;
    cout << Q.rows() << " " << Q.cols() << endl;
    //cout << Q.transpose() * Q << endl;
    cout << "it worked" << endl << endl;
    //Eigen::MatrixXd Hp = Q.transpose() * (H + 0.0 * MatrixXd::Identity(dim, dim)) * Q;
    //Eigen::MatrixXd Sp = Q.transpose() * (S + 0.0 * MatrixXd::Identity(dim, dim)) * Q;
    Eigen::MatrixXd Hp = Q.transpose() * H * Q;
    Eigen::MatrixXd Sp = Q.transpose() * S * Q;
    //cout << Hp << endl << endl;
    //cout << Sp << endl << endl;
    Eigen::GeneralizedEigenSolver<MatrixXd> es1(Hp, Sp);
    cout << "Lowest eigenpair of H, S" << endl << endl;;
    VectorXcd Lambda = es1.eigenvalues();
    MatrixXcd V = Q * es1.eigenvectors();
    SortEig(Lambda, V);
    std::complex<double> val = Lambda(0);
    Eigen::VectorXcd v = V.col(0);
    cout << val << endl << endl;
    //cout << v.transpose() << endl << endl;
    cout << "residual" << endl;
    cout << ((H - val * S) * v).norm() << endl << endl;
    cout << "all eigenvalues" << endl;
    for (int i = 0; i < 10; i++) { cout << Lambda(i) << " "; }
    cout << endl << endl;
    //cout << Lambda.transpose() << endl << endl;

    auto endRand = std::chrono::steady_clock::now();
    std::cout << "Time for RandSolver: " << std::chrono::duration_cast<std::chrono::seconds>(endRand - beginRand).count() << "[s]" << std::endl;

    //std::ofstream f("vec" + std::string(argv[1]) + ".txt");
    std::ofstream f("vec.txt");
    for (int i = 0; i < v.size(); i++)
    {
        //f << lowvec(i) / lowvec(0) << " | " << v(i) / v(0) << endl;
        f << lowvec(i) << " | " << v(i) << endl;
    }

}
