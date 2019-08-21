#ifndef DAVIDSON_HEADER
#define DAVIDSON_HEADER
#include <numeric>
#include <Eigen/Dense>

void OrthonormalizeVectorToSubspace(Eigen::VectorXd &delta, const Eigen::MatrixXd &V)
{
    //modified grahm schmidt
    delta.normalize();
    for (int i = 0; i < V.cols(); i++)
    {
        double alpha = V.col(i).dot(delta);
        delta = delta - V.col(i) * alpha;
        delta.normalize();
    }
}

void SortEig(Eigen::VectorXd &D, Eigen::MatrixXd &V)
{
    Eigen::VectorXd Dcopy(D);
    Eigen::MatrixXd Vcopy(V);
    std::vector<double> a(Dcopy.data(), Dcopy.data() + Dcopy.size());
    //initialize original index locations
    std::vector<int> idx(a.size());
    std::iota(idx.begin(), idx.end(), 0);
    //sort indexes based on comparing values in a
    std::sort(idx.begin(), idx.end(), [&a](int i, int j) -> bool { return a[i] < a[j]; });
    for (int i = 0; i < idx.size(); i++)
    {
        D(i) = Dcopy(idx[i]);
        V.col(i) = Vcopy.col(idx[i]);
    }
}

int FindClosest(double val, const Eigen::VectorXd &V)
{
    std::vector<double> a(V.data(), V.data() + V.size());
    auto begin = a.begin();
    auto end = a.end();
    auto it = std::upper_bound(begin, end, val);
    if (it == begin) { return 0; }
    if (it == end) { return std::distance(begin, end) - 1; }
    double diff1 = std::abs(*it - val);
    double diff2 = std::abs(*(it - 1) - val);
    if (diff2 < diff1) { --it; }
    return std::distance(begin, it);
}

class Davidson
{
    protected:
    int n, vmax, nrestart;
    Eigen::MatrixXd Vectors;
    Eigen::VectorXd Values;

    public:
    Davidson() {}
    Davidson(const Eigen::MatrixXd &A, int _n = 1, int _vmax = 25, int _nrestart = 5) : n{_n}, vmax{_vmax}, nrestart{_nrestart} { run(A, n, vmax, nrestart); }

    Eigen::VectorXd eigenvalues() const { return Values; }
    Eigen::MatrixXd eigenvectors() const { return Vectors; }

    int run(const Eigen::MatrixXd &A, int _n = 1, int _vmax = 25, int _nrestart = 5) 
    {
        n = _n;
        vmax = _vmax;
        nrestart = _nrestart;
        double tol = 1.e-8;
        int dim = A.rows();
        assert(A.rows() == A.cols());
        //if small problem
        if (dim < vmax)
        {
            vmax = dim;
            int a = 0.1 * dim;
            nrestart = a > 1 ? a : 1;
        }
        Vectors.setZero(dim, n);
        Values.setZero(n);

        //initalize vector and subspace
        int i;
        A.diagonal().minCoeff(&i);
        Eigen::VectorXd z = Eigen::VectorXd::Unit(dim, i);
        Eigen::VectorXd Az = A * z;
        Eigen::MatrixXd V(z), AV(Az);

        double target = z.adjoint() * Az;
        double old_theta = 1000000;
        double old_rNorm = 1000000;
        for (int iter = 1; iter <= 1000; iter++)
        {    
            //subspace problem
            Eigen::MatrixXd Aprime = V.adjoint() * AV;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Aprime);
            Eigen::VectorXd D = es.eigenvalues();
            Eigen::MatrixXd U = es.eigenvectors();

            //sort solution and extract answer
            SortEig(D, U);
            int index = FindClosest(target, D);
            double theta = D(index);
        std::cout << theta << std::endl;
            Eigen::VectorXd s = U.col(index);
            Eigen::VectorXd u = V * s;
            Eigen::VectorXd Au = AV * s;

            //residual
            Eigen::VectorXd r = Au - theta * u;
            double rNorm = r.norm();

        std::cout << std::abs(theta - old_theta) << std::endl;
            //check for convergence
            if (std::abs(theta - old_theta) < tol)
            {
                Values(0) = theta;
                Vectors.col(0) = u;
                return iter;
            }
            old_theta = theta;

            //restart
            if (V.cols() >= vmax)
            {
                Eigen::MatrixXd N = U.block(0, index, U.rows(), nrestart);
                V = V * N;
                AV = AV * N;
            }

            //update target
            if (rNorm < old_rNorm)
            {
                target = theta;
                old_rNorm = rNorm;
            }

            //calculate correction
            Eigen::VectorXd delta(r.rows());
            for (int l = 0; l < delta.rows(); l++)
            {
                delta(l) = r(l) / (A(l, l) - theta + 1.e-8);
            }
            OrthonormalizeVectorToSubspace(delta, V);

            //append vector to subspace
            int m = V.cols();
            V.conservativeResize(dim, m + 1);
            AV.conservativeResize(dim, m + 1);
            V.col(m) = delta;
            AV.col(m) = A * delta;
        }
        //failed
        return -1;
    }

    int run(const Operator::Hamiltonian &A, const FockVector &CI, int _n = 1, int _vmax = 25, int _nrestart = 5) 
    {
        n = _n;
        vmax = _vmax;
        nrestart = _nrestart;
        double tol = 1.e-8;
        int dim = CI.size();
        //diagonal
        Eigen::VectorXd diag;
        A.diagonal(CI, diag);
        //if small problem
        if (dim < vmax)
        {
            vmax = dim;
            int a = 0.1 * dim;
            nrestart = a > 1 ? a : 1;
        }
        Vectors.setZero(dim, n);
        Values.setZero(n);

        //initalize vector and subspace
        int i;
        diag.minCoeff(&i);
        Eigen::VectorXd z = Eigen::VectorXd::Unit(dim, i);
        //Eigen::VectorXd Az = A * z;
        Eigen::VectorXd Az = A.multiply(CI, z);
        Eigen::MatrixXd V(z), AV(Az);

        double target = z.adjoint() * Az;
        double old_theta = 1000000;
        double old_rNorm = 1000000;
        for (int iter = 1; iter <= 1000; iter++)
        {    
            //subspace problem
            Eigen::MatrixXd Aprime = V.adjoint() * AV;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Aprime);
            Eigen::VectorXd D = es.eigenvalues();
            Eigen::MatrixXd U = es.eigenvectors();

            //sort solution and extract answer
            SortEig(D, U);
            int index = FindClosest(target, D);
            double theta = D(index);
        std::cout << theta << std::endl;
            Eigen::VectorXd s = U.col(index);
            Eigen::VectorXd u = V * s;
            Eigen::VectorXd Au = AV * s;

            //residual
            Eigen::VectorXd r = Au - theta * u;
            double rNorm = r.norm();

            //check for convergence
        std::cout << std::abs(theta - old_theta) << std::endl;
            if (std::abs(theta - old_theta) < tol)
            {
                Values(0) = theta;
                Vectors.col(0) = u;
                return iter;
            }
            old_theta = theta;

            //restart
            if (V.cols() >= vmax)
            {
                Eigen::MatrixXd N = U.block(0, index, U.rows(), nrestart);
                V = V * N;
                AV = AV * N;
            }

            //update target
            if (rNorm < old_rNorm)
            {
                target = theta;
                old_rNorm = rNorm;
            }

            //calculate correction
            Eigen::VectorXd delta(r.rows());
            for (int l = 0; l < delta.rows(); l++)
            {
                delta(l) = r(l) / (diag(l) - theta + 1.e-8);
            }
            OrthonormalizeVectorToSubspace(delta, V);

            //append vector to subspace
            int m = V.cols();
            V.conservativeResize(dim, m + 1);
            AV.conservativeResize(dim, m + 1);
            V.col(m) = delta;
            //AV.col(m) = A * delta;
            AV.col(m) = A.multiply(CI, delta);
        }
        //failed
        return -1;
    }
};
#endif
