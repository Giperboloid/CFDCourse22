
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
        InitMatrixWithoutOwning(A, SystemStencil);

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

        InitMatrixWithoutOwning(D, DiagonalStencil);

        GetTriangularComponent(UpperStencil, SystemStencil);
        GetTriangularComponent(LowerStencil, SystemStencil, false);

        InitMatrixWithoutOwning(U, UpperStencil);
        InitMatrixWithoutOwning(L, LowerStencil);

        SetOptionalRelaxParameter();
        if(RelaxParam != 1)
            throw std::runtime_error("START TEST");

    }

    /**
     *
     * @return
     */
    bool CheckSystemForSymmetry()
    {
        return false;
    }

    /**
     *
     */
    void SolveSor() {};

    /**
     *
     */
    void SolveSSor() {};

private:

    /**
    *
    * @param matrix
    * @param stencil_ptr
    */
    static void InitMatrixWithoutOwning(Matrix& matrix, Stencil& stencil_ptr) {
        matrix.own_data = false;
        matrix.nrows = matrix.ncols = stencil_ptr.n_rows;
        matrix.nnz = stencil_ptr.nnz;
        matrix.ptr = stencil_ptr.addr.data();
        matrix.col = stencil_ptr.cols.data();
        matrix.val = stencil_ptr.vals.data();
    }

    /**
     *
     * @param out_stencil
     * @param init_stencil
     * @param upper
     */
    static void GetTriangularComponent(Stencil& out_stencil, const Stencil& init_stencil, bool upper = true) {

        /// Creating upper triangular matrix U:
        std::vector<long> U_addr(init_stencil.addr.size());
        std::vector<long> U_cols;
        std::vector<double> U_vals;

        // just fill with old (system matrices) values here
        for(auto i = 0; i < init_stencil.addr.size(); ++i)
            U_addr[i] = init_stencil.addr[i];

        for(auto i = 0; i < init_stencil.addr.size() - 1; ++i)
        {
            for(auto j = init_stencil.addr[i]; j < init_stencil.addr[i + 1]; ++j)
            {
                if(upper)
                {
                    if(init_stencil.cols[j] > i)
                    {
                        U_vals.emplace_back(init_stencil.vals[j]);
                        U_cols.emplace_back(init_stencil.cols[j]);
                    }
                    else
                        for(auto k = i; k < U_addr.size() - 1; ++k, U_addr[k]-=1);
                }
                else
                {
                    if(init_stencil.cols[j] < i)
                    {
                        U_vals.emplace_back(init_stencil.vals[j]);
                        U_cols.emplace_back(init_stencil.cols[j]);
                    }
                    else
                        for(auto k = i; k < U_addr.size() - 1; ++k, U_addr[k]-=1);
                }
            }
        }

        amgcl::backend::numa_vector<long> numa_U_addr (U_addr);
        amgcl::backend::numa_vector<long> numa_U_cols (U_cols);
        amgcl::backend::numa_vector<double> numa_U_vals (U_vals);

        out_stencil.addr.swap(numa_U_addr);
        out_stencil.cols.swap(numa_U_cols);
        out_stencil.vals.swap(numa_U_vals);

        out_stencil.n_rows = init_stencil.n_rows;
        out_stencil.nnz = static_cast<long>(out_stencil.vals.size());
    }

    /**
     *
     */
    void SetOptionalRelaxParameter() {

        auto triangles_sum_ptr = amgcl::backend::sum(1.0, U, 1.0, L);
        if(!triangles_sum_ptr)
            throw std::runtime_error("Setting relaxation parameter error: can not get sum of triangle matrices");

        amgcl::backend::numa_vector<double>Inverse_D_vals(DiagonalStencil.n_rows);

        for(auto i = 0; i < Inverse_D_vals.size(); ++i)
            Inverse_D_vals[i] = 1.0 / DiagonalStencil.vals[i];

        /// Temporary set D matrices values pointer to 'Inverse_D_vals' <-> gives opportunity to use D as D^-1
        /// without creating new one matrix
        D.val = Inverse_D_vals.data();

        auto product_matrix = amgcl::backend::product(D, *triangles_sum_ptr);
        if(!product_matrix)
            throw std::runtime_error("Setting relaxation parameter error: can not get matrices product");

        double spectral_radius = amgcl::backend::spectral_radius<false, Matrix>(*product_matrix);

        RelaxParam = 2.0 / (1.0 + sqrt(1.0 - spectral_radius * spectral_radius));

        /// Revert D matrices values
        D.val = DiagonalStencil.vals.data();
    }


    Stencil SystemStencil, DiagonalStencil,
            LowerStencil, UpperStencil;

    /// A = D + L + U factorization parts:
    Matrix D, L, U;

    /// Optional relaxation parameter: (0, 2); initially 1.
    double RelaxParam {1};


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