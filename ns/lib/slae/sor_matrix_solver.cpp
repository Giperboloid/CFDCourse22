
#include "sor_matrix_solver.hpp"
#include <amgcl/backend/builtin.hpp>

struct SorMatrixSolver::Implementation {

public:

    using Matrix = amgcl::backend::crs<double, int>;

    explicit Implementation(const Matrix& matrix)
    : Diagonal{amgcl::backend::diagonal(matrix)}
    {
        if(!CheckConvergenceCondition())
            throw std::runtime_error("Necessary condition is not met: solution will not converge");
    }

    bool CheckConvergenceCondition()
    {

    }

    bool CheckForSymmetry()
    {

    }

    void SolveSor() {};
    void SolveSSor() {};

private:

    /// LU factorization elements
    std::shared_ptr<amgcl::backend::numa_vector<double>> Diagonal;
    Matrix L, U;

    uint8_t RelaxParam {1};


};

SorMatrixSolver::SorMatrixSolver(unsigned MaxIter, double Tol, bool IsSym)
: MaxIterations(MaxIter), Tolerance(Tol), IsInitiallySymmetric(IsSym)
{}

void SorMatrixSolver::set_matrix(const CsrStencil &stencil, const std::vector<double> &mat_values)
{
    Implementation::Matrix matrix;
    matrix.own_data = false;
    matrix.nrows = matrix.ncols = stencil.n_rows();
    matrix.nnz = stencil.n_nonzero();
    matrix.ptr = const_cast<int*>(stencil.addr().data());
    matrix.col = const_cast<int*>(stencil.cols().data());
    matrix.val = const_cast<double*>(mat_values.data());


    Impl = std::make_unique<SorMatrixSolver::Implementation>(matrix);

}

void SorMatrixSolver::solve(const std::vector<double> &rhs, std::vector<double> &ret) const {
    if(Impl)
    {
        if(IsInitiallySymmetric || Impl->CheckForSymmetry())
            Impl->SolveSSor();
        else
            Impl->SolveSor();
    }
    else
        throw std::runtime_error("Matrix was not passed to the solver");
}