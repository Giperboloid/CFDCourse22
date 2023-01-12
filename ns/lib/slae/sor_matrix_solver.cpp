
#include "sor_matrix_solver.hpp"

#include <amgcl/backend/builtin.hpp>

class SorMatrixSolver::Implementation {

    using Matrix = amgcl::backend::crs<double, long>;

public:

    typedef struct MatrixStencil {
        amgcl::backend::numa_vector<long> addr;
        amgcl::backend::numa_vector<long> cols;
        amgcl::backend::numa_vector<double> vals;
        long n_rows{0};
        long nnz{0};
    } Stencil;

    Implementation(const CsrStencil &stencil, const std::vector<double> &mat_values)
    : SystemStencil {
        .addr = std::vector<long>(stencil.addr().begin(), stencil.addr().end()),
        .cols = std::vector<long>(stencil.cols().begin(), stencil.cols().end()),
        .vals = mat_values,
        .n_rows = stencil.n_rows(),
        .nnz = stencil.n_nonzero()
    }
    {}

    /**
     * @brief Following method implements A = D + L + U factorization procedure:
     * stencils creation and matrices filling
     */
    void Initialize() {

        /// Set locally system matrix A:
        Matrix A;
        InitMatrixWithoutOwning(A, &SystemStencil);

        /// Set D matrix - system matrices diagonal component:
        auto D_vals = amgcl::backend::diagonal(A);
        if(!D_vals)
            throw std::runtime_error("Solver init error: cant get system's diagonal");

        (*D_vals).swap(DiagonalStencil.vals);

        if(DiagonalStencil.vals.size() != A.nrows)
            throw std::runtime_error("Solver init error: the main diagonal of system matrix contains zeroes, "
                                     "SOR(Symmetric SOR) method can not be used");

        DiagonalStencil.n_rows = DiagonalStencil.nnz = static_cast<long>(A.nrows);
        DiagonalStencil.addr.resize(DiagonalStencil.n_rows + 1);
        DiagonalStencil.cols.resize(DiagonalStencil.n_rows);


        for(long i = 0; i <= DiagonalStencil.n_rows; ++i) {
            if(i != DiagonalStencil.n_rows)
            {
                DiagonalStencil.addr[i] = i;
                DiagonalStencil.cols[i] = i;
            }
            else
                DiagonalStencil.addr[i] = i;
        }

        InitMatrixWithoutOwning(D, &DiagonalStencil);

        SetOptionalRelaxParameter();
        if(!RelaxParam)
            throw std::runtime_error("START TEST");

    }

    void SetOptionalRelaxParameter() {
        RelaxParam = 0;
    }

    bool CheckSystemForSymmetry()
    {
        return false;
    }

    static void InitMatrixWithoutOwning(Matrix& matrix, Stencil* stencil_ptr) {
        matrix.own_data = false;
        matrix.nrows = matrix.ncols = stencil_ptr->n_rows;
        matrix.nnz = stencil_ptr->nnz;
        matrix.ptr = stencil_ptr->addr.data();
        matrix.col = stencil_ptr->cols.data();
        matrix.val = stencil_ptr->vals.data();
    }

    void SolveSor() {};
    void SolveSSor() {};

private:

    Stencil SystemStencil, DiagonalStencil,
            LowerStencil, UpperStencil;

    /// A = D + L + U factorization parts:
    Matrix D, L, U;

    /// Optional relaxation parameter: (0, 2); initially 1.
    uint8_t RelaxParam {1};


};

SorMatrixSolver::SorMatrixSolver(unsigned MaxIter, double Tol, bool IsSym)
: MaxIterations(MaxIter), Tolerance(Tol), IsInitiallySymmetric(IsSym)
{}

void SorMatrixSolver::set_matrix(const CsrStencil &stencil, const std::vector<double> &mat_values)
{
    Impl = std::make_unique<SorMatrixSolver::Implementation>(stencil, mat_values);
    Impl->Initialize();
}

void SorMatrixSolver::solve(const std::vector<double> &rhs, std::vector<double> &ret) const {
    if(Impl)
    {
        if(IsInitiallySymmetric || Impl->CheckSystemForSymmetry())
            Impl->SolveSSor();
        else
            Impl->SolveSor();
    }
    else
        throw std::runtime_error("Matrix was not passed to the solver");
}