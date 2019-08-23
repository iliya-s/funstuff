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

//returns index of vector element that is closest to val
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

//returns index of vector element that is just lower than val
int FindLower(double val, const Eigen::VectorXd &V)
{
    std::vector<double> a(V.data(), V.data() + V.size());
    auto begin = a.begin();
    auto end = a.end();
    auto it = std::upper_bound(begin, end, val);
    if (it == begin) { return 0; }
    return std::distance(begin, it - 1);
} 

//abstract class for multiplying a matrix onto a vector, makes it so only one davidson implementation is needed
class AbstractMatrixMult
{
    private:
    public:
    virtual int dimension() const = 0;
    virtual void diagonal(Eigen::VectorXd &V) const = 0;
    virtual void multiply(const Eigen::VectorXd &V, Eigen::VectorXd &AV) const = 0;
};

class MatrixMult : public AbstractMatrixMult
{
    private:
    Eigen::MatrixXd H;

    public:
    MatrixMult(const Operator::Hamiltonian &Ham, const FockVector &V) { Ham.matrix(V, H); }
    MatrixMult(const Eigen::MatrixXd &Ham) : H{Ham} { assert(Ham.rows() == Ham.cols()); }
    int dimension() const { return H.rows(); }
    void diagonal(Eigen::VectorXd &V) const { V = H.diagonal(); }
    void multiply(const Eigen::VectorXd &V, Eigen::VectorXd &HV) const { HV = H * V; }
};

class DirectMatrixMult : public AbstractMatrixMult
{
    private:
    const Operator::Hamiltonian &H;
    const FockVector &CI;

    public:
    DirectMatrixMult(const Operator::Hamiltonian &Ham, const FockVector &V) : H{Ham}, CI{V} {}
    int dimension() const { return CI.size(); }
    void diagonal(Eigen::VectorXd &V) const { H.diagonal(CI, V); }
    void multiply(const Eigen::VectorXd &V, Eigen::VectorXd &HV) const { H.multiply(CI, V, HV); }
};

class Davidson
{
    private:
    int n, vmax, nrestart;
    double tol;
    Eigen::MatrixXd Vectors;
    Eigen::VectorXd Values;

    public:
    Davidson() {}
    Davidson(const Eigen::MatrixXd &A, int _n = 1, int _vmax = 25, int _nrestart = 5, double _tol = 1.e-8) : n{_n}, vmax{_vmax}, nrestart{_nrestart}, tol{_tol} { run(A, n, vmax, nrestart, tol); }
    Davidson(const AbstractMatrixMult &A, int _n = 1, int _vmax = 25, int _nrestart = 5, double _tol = 1.e-8) : n{_n}, vmax{_vmax}, nrestart{_nrestart}, tol{_tol} { run(A, n, vmax, nrestart, tol); }

    const Eigen::VectorXd &eigenvalues() const { return Values; }
    const Eigen::MatrixXd &eigenvectors() const { return Vectors; }

    int run(const Eigen::MatrixXd &A, int _n = 1, int _vmax = 25, int _nrestart = 5, double _tol = 1.e-8) 
    {
        MatrixMult Mat(A);
        return run(Mat, _n, _vmax, _nrestart, _tol);
    }

    int run(const Operator::Hamiltonian &H, const FockVector &CI, int _n = 1, int _vmax = 25, int _nrestart = 5, double _tol = 1.e-8) 
    {
        DirectMatrixMult Mat(H, CI);
        return run(Mat, _n, _vmax, _nrestart, _tol);
    }

    int run(const AbstractMatrixMult &A, int _n = 1, int _vmax = 25, int _nrestart = 5, double _tol = 1.e-8) 
    {
        n = _n;
        vmax = _vmax;
        nrestart = _nrestart;
        tol = _tol;
        int dim = A.dimension();
        //diagonal
        Eigen::VectorXd diag;
        A.diagonal(diag);
        //if small problem
        if (dim < vmax)
        {
            vmax = dim;
            int a = 0.1 * dim;
            nrestart = a > 1 ? a : 1;
        }
        //make room for solution
        Vectors.setZero(dim, n);
        Values.setZero(n);

        //initalize vector and subspace
        int i;
        diag.minCoeff(&i);
        Eigen::VectorXd z = Eigen::VectorXd::Unit(dim, i);
        Eigen::VectorXd Az;
        A.multiply(z, Az);
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
            //int index = FindClosest(target, D);
            int index = FindLower(target, D);
            double theta = D(index);
            Eigen::VectorXd s = U.col(index);
            Eigen::VectorXd u = V * s;
            Eigen::VectorXd Au = AV * s;

            //residual
            Eigen::VectorXd r = Au - theta * u;
            double rNorm = r.norm();

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
                delta(l) = r(l) / (diag(l) - theta + 1.e-8);
            }
            OrthonormalizeVectorToSubspace(delta, V);

            //append vector to subspace
            int m = V.cols();
            V.conservativeResize(dim, m + 1);
            AV.conservativeResize(dim, m + 1);
            V.col(m) = delta;
            Eigen::VectorXd Hdelta;
            A.multiply(delta, Hdelta);
            AV.col(m) = Hdelta;
        }
        //failed
        return -1;
    }
};
#endif
