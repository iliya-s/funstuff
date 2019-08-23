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
            assert(n < dim);
        }
        //make room for solution
        Vectors.setZero(dim, n);
        Values.setZero(n);

        //initalize vector and subspace
        Eigen::MatrixXd V(dim, n), AV(dim, n);
        Eigen::VectorXd diagcopy(diag);
        for (int i = 0; i < n; i++)
        {
            int index;
            diagcopy.minCoeff(&index);
            Eigen::VectorXd z = Eigen::VectorXd::Unit(dim, index);
            Eigen::VectorXd Az = Eigen::VectorXd::Zero(dim);
            //Az(index) = diagcopy(index);
            A.multiply(z, Az);
            V.col(i) = z;
            AV.col(i) = Az;
            diagcopy(index) = 1000000;
        }
        double target = V.col(0).adjoint() * AV.col(0); //target will always be ground state

        Eigen::VectorXd old_Theta = Eigen::VectorXd::Constant(n, 1000000);
        double old_rNorm = 1000000;
        for (int iter = 1; iter <= 1000; iter++)
        {    
            //std::cout << iter << std::endl;
            //subspace problem
            int subdim = V.cols(); //dimension of subspace
            Eigen::MatrixXd Aprime = V.adjoint() * AV;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Aprime);
            Eigen::VectorXd D = es.eigenvalues();
            Eigen::MatrixXd S = es.eigenvectors();

            //sort solution and extract answer
            SortEig(D, S);
            int index = FindClosest(target, D);
            //int index = FindLower(target, D);
            //int index = 0;
            
            Eigen::VectorXd Theta = D.segment(index, n);
            Eigen::MatrixXd U = V * S.block(0, index, subdim, n);
            Eigen::MatrixXd AU = AV * S.block(0, index, subdim, n);
            Eigen::MatrixXd Delta = Eigen::MatrixXd::Zero(dim, n);
            Eigen::MatrixXd HDelta = Eigen::MatrixXd::Zero(dim, n);

            for (int i = 0; i < n; i++)
            {
                //residual
                Eigen::VectorXd r = AU.col(i) - Theta(i) * U.col(i);

                if (i == 0) //target will always be ground state
                {
                    //update target
                    double rNorm = r.norm();
                    if (rNorm < old_rNorm)
                    {
                        target = Theta(i);
                        old_rNorm = rNorm;
                    }
                }
                
                //calculate correction
                Eigen::VectorXd delta(dim);
                for (int l = 0; l < delta.rows(); l++) { delta(l) = - r(l) / (diag(l) - Theta(i) + 1.e-8); }
                OrthonormalizeVectorToSubspace(delta, V);
                Eigen::VectorXd Hdelta;
                A.multiply(delta, Hdelta);
                Delta.col(i) = delta;
                HDelta.col(i) = Hdelta;
            }

            //grow subspace and append vectors
            V.conservativeResize(dim, subdim + n);
            AV.conservativeResize(dim, subdim + n);
            V.block(0, subdim, dim, n) = Delta;
            AV.block(0, subdim, dim, n) = HDelta;

            //check for convergence
            if ((Theta - old_Theta).squaredNorm() < tol)
            {
                Values = Theta;
                Vectors = U;
                return iter;
            }
            old_Theta = Theta;

            //restart
            if (V.cols() >= vmax)
            {
                Eigen::MatrixXd N = U.block(0, index, dim, nrestart);
                V = V * N;
                AV = AV * N;
            }
        }
        //failed
        return -1;
    }
};
#endif
