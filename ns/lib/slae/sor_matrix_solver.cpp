
#include "sor_matrix_solver.hpp"

#include <amgcl/backend/builtin.hpp>

struct SorMatrixSolver::Implementation {

public:

    using Matrix = amgcl::backend::crs<double, long>;

    Implementation(const CsrStencil &stencil, const std::vector<double> &mat_values)
    {
        /// Set system matrix A:
        A.own_data = false;
        A.nrows = A.ncols = stencil.n_rows();
        A.nnz = stencil.n_nonzero();

        std::vector<long> stencil_addr_long(stencil.addr().begin(), stencil.addr().end());
        std::vector<long> stencil_col_long(stencil.cols().begin(), stencil.cols().end());

        A.ptr = const_cast<long*>(stencil_addr_long.data());
        A.col = const_cast<long*>(stencil_col_long.data());

        A.val = const_cast<double*>(mat_values.data());

        /// Set D matrix - system matrices diagonal component:
        auto D_vector = amgcl::backend::diagonal(A);

        if(D_vector->size() != A.nrows)
            throw std::runtime_error("The main diagonal of system matrix contains zeroes: "
                                     "SOR(SSOR) method can not be used");

        D.own_data = false;
        D.nrows = D.ncols = D.nnz = stencil.n_rows();
        D.val = const_cast<double *>(D_vector->data());

        std::vector<long> D_addr_vector{};
        std::vector<long> D_col_vector{};

        for(long i = 0; i <= static_cast<long>(D.ncols); ++i) {
            if(i != D.ncols) {
                D_addr_vector.emplace_back(i);
                D_col_vector.emplace_back(i);
            }
            else
                D_addr_vector.emplace_back(i);
        }
        D.ptr = const_cast<long *>(D_addr_vector.data());
        D.col = const_cast<long *>(D_col_vector.data());

        /// Set D^-1 matrix - inverse system matrices diagonal component:
        auto D_inverse_vector = amgcl::backend::diagonal(A, true);
        Inverse_D = D;
        Inverse_D.val = const_cast<double *>(D_inverse_vector->data());


        if(!CheckConvergenceCondition())
            throw std::runtime_error("Necessary condition is not met: solution will not converge");


    }

    bool CheckConvergenceCondition()
    {
        /// Prepare a unit matrix:
        Matrix I {D};
        std::vector<double> I_vector (I.nrows, 1);
        I.val = const_cast<double *>(I_vector.data());

        /// Check necessary convergence condition: Mu = p(C) = p(I - (D^-1)A)
        std::shared_ptr<Matrix> Inverse_D_product_A = amgcl::backend::product(Inverse_D, A);

        std::shared_ptr<Matrix> Jacobi_Iter_Matrix = amgcl::backend::sum(1.0, I, -1.0, *Inverse_D_product_A);

        Mu = amgcl::backend::spectral_radius<false, Matrix>(*Jacobi_Iter_Matrix);

        if(Mu < 1)
            return true;

        return false;
    }

    bool CheckForSymmetry()
    {
        return false;
    }

    void SolveSor() {};
    void SolveSSor() {};

private:

    /// LU factorization parts
    Matrix A, D, Inverse_D, L, U;

    /// Spectral radius of Jacobi's iteration matrix C = I - (D^-1)A
    double Mu {0};

    /// Optional relaxation parameter: (0, 2)
    uint8_t RelaxParam {1};


};

SorMatrixSolver::SorMatrixSolver(unsigned MaxIter, double Tol, bool IsSym)
: MaxIterations(MaxIter), Tolerance(Tol), IsInitiallySymmetric(IsSym)
{}

void SorMatrixSolver::set_matrix(const CsrStencil &stencil, const std::vector<double> &mat_values)
{
    Impl = std::make_unique<SorMatrixSolver::Implementation>(stencil, mat_values);
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