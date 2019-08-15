#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

using namespace Eigen;
using namespace std;

void GeneralizedJacobiDavidson(const Eigen::MatrixXd &H, const Eigen::MatrixXd &S, double &lambda, Eigen::VectorXd &v)
{
    int dim = H.rows(); //dimension of problem
    int restart = std::max((int) (0.1 * dim), 20); //20 - 30
    int q = std::max((int) (0.04 * dim), 5); //5 - 10
    Eigen::MatrixXd V, HV, SV;  //matrices storing action of vector on sample space
    //Eigen::VectorXd z = Eigen::VectorXd::Random(dim);
    Eigen::VectorXd z = Eigen::VectorXd::Unit(dim, 0);
    while (1)
    {
        int m = V.cols(); //number of vectors in subspace at current iteration

        //modified grahm schmidt to orthogonalize sample space with respect to overlap matrix
        Eigen::VectorXd Sz = S * z;
        for (int i = 0; i < m; i++)
        {
            double alpha = V.col(i).adjoint() * Sz;
            z = z - alpha * V.col(i);
        }

        //normalize z after orthogonalization and calculate action of matrices on z
        Sz = S * z;
        double beta = std::sqrt(z.adjoint() * Sz);
        Eigen::VectorXd v_m = z / beta;
        Eigen::VectorXd Hv_m = H * v_m;
        Eigen::VectorXd Sv_m = S * v_m;
        //store new vector
        V.conservativeResize(H.rows(), m + 1);
        HV.conservativeResize(H.rows(), m + 1);
        SV.conservativeResize(H.rows(), m + 1);
        V.col(m) = v_m;
        HV.col(m) = Hv_m;
        SV.col(m) = Sv_m;

        //solve eigenproblem in subspace
        Eigen::MatrixXd A = V.adjoint() * HV;
        Eigen::MatrixXd B = V.adjoint() * SV;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, B);
        //if sample space size is too large, restart subspace
        if (m == restart)
        {
            Eigen::MatrixXd N = es.eigenvectors().block(0, 0, A.rows(), q);
            V = V * N;
            HV = HV * N;
            SV = SV * N;
            A = V.adjoint() * HV;
            B = V.adjoint() * SV;
            es.compute(A, B);
        }
        Eigen::VectorXd s = es.eigenvectors().col(0);
        double theta = es.eigenvalues()(0);

        cout << theta << "  " << m << endl;
        //transform vector into original space
        Eigen::VectorXd u = V * s;
        //calculate residue vector
        Eigen::VectorXd u_H = HV * s;
        Eigen::VectorXd u_S = SV * s;
        Eigen::VectorXd r = u_H - theta * u_S;

        if (r.squaredNorm() < 1.e-8)
        {
            lambda = theta;
            v = u;
            break;
        }
        
        Eigen::MatrixXd X = Eigen::MatrixXd::Identity(u.rows(), u.rows());
        Eigen::MatrixXd Xleft = X - S * u * u.adjoint();
        Eigen::MatrixXd Xright = X - u * u.adjoint() * S;
        X = Xleft * (H - theta * S) * Xright;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(X);
        z = dec.solve(-r);
    }
}

int main(int argc, char **argv)
{
    /*
    int D = std::atoi(argv[1]);
    MatrixXd X = MatrixXd::Random(D, D);
    MatrixXd H = X + X.adjoint();

    X = MatrixXd::Random(D, D);
    MatrixXd S = X + X.adjoint();
    S += D * MatrixXd::Identity(D, D);
    */
    MatrixXd H, S;
    ifstream Hin("H.txt");
    ifstream Sin("S.txt");
    int D = std::count(std::istreambuf_iterator<char>(Hin), std::istreambuf_iterator<char>(), '\n');
    Hin.seekg(0, ios::beg);
    H.resize(D, D);
    S.resize(D, D);
    for (int i = 0; i < D; i++)
    {
        for (int j = 0; j < D; j++)
        {
            Hin >> H(i, j);
            Sin >> S(i, j);
        }
    }
    Hin.close();
    Sin.close();

    cout << "H, S" << endl << endl;
    cout << H << endl << endl;
    cout << S << endl << endl;

    MatrixXd Sprime = S + 0.001 * MatrixXd::Identity(D, D);
    GeneralizedSelfAdjointEigenSolver<MatrixXd> es(H, Sprime);

    cout << "Lowest eigenpair of H, S" << endl << endl;
    VectorXd lowvec = es.eigenvectors().col(0);
    cout << es.eigenvalues()(0) << endl << endl;
    //cout << es.eigenvectors().col(0) << endl << endl;
    
    cout << "Using Davidson" << endl << endl;
    double theta;
    Eigen::VectorXd v;
    GeneralizedJacobiDavidson(H, S, theta, v);

    cout << theta << endl << endl;
    for (int i = 0; i < lowvec.size(); i++)
    {
        cout << lowvec(i) << "    " << v(i) << endl;
    }
    cout << endl << endl << lowvec.adjoint() * S * lowvec << " " << v.adjoint() * S * v << endl;
    cout << endl << lowvec.adjoint() * H * lowvec << " " << v.adjoint() * H * v << endl;
    
    /*`
    for (int i = 1; i < lowvec.size(); i++)
    {
        cout << lowvec(i)/lowvec(0) << "    " << v(i)/v(0) << endl;
    }
    */

    /*
    cout << "trying to block it" << endl;

    MatrixXd H1, H2, S1, S2;
    //top left block
    H1 = H.block(0, 0, D / 2, D /2);
    S1 = S.block(0, 0, D / 2, D / 2);

    //bottom right block
    H2 = H.block(D / 2, D /2, D / 2, D / 2);
    S2 = S.block(D / 2, D / 2, D / 2, D / 2);

    cout << H1 << endl << endl;
    cout << H2 << endl << endl;

    es.compute(H1, S1);
    VectorXd d1 = es.eigenvalues();
    MatrixXd v1 = es.eigenvectors();

    cout << es.eigenvalues() << endl << endl;
    cout << es.eigenvectors() << endl << endl;

    es.compute(H2, S2);
    VectorXd d2 = es.eigenvalues();
    MatrixXd v2 = es.eigenvectors();

    cout << es.eigenvalues() << endl << endl;
    cout << es.eigenvectors() << endl << endl;

    MatrixXd V(D, 2);
    V.col(0) << v1.col(0), v2.col(0);
    V.col(1) << v1.col(1), v2.col(1);

    MatrixXd H_prime = V.adjoint() * H * V;
    MatrixXd S_prime = V.adjoint() * S * V;
    cout << "subspace problem" << endl;
    cout << H_prime << endl << endl;
    cout << S_prime << endl << endl;
    es.compute(H_prime, S_prime);

    cout << es.eigenvalues() << endl << endl;
    cout << es.eigenvectors() << endl << endl;

    VectorXd answer = V * es.eigenvectors().col(0);

    for (int i = 1; i < lowvec.size() - 1; i++)
    {
        cout << lowvec(i)/lowvec(0) << "\t" << answer(i)/answer(0) << "\t" << lowvec(i)/lowvec(0) - answer(i)/answer(0) << endl;
    }
    */
}
