
#include "sor_matrix_solver.hpp"

#include <amgcl/backend/builtin.hpp>

class SorMatrixSolver::Implementation {

    using Matrix = amgcl::backend::crs<double, long>;

    const double UniversalRelaxParameter = 1.95;

    typedef struct MatrixStencil {
        amgcl::backend::numa_vector<long> addr;
        amgcl::backend::numa_vector<long> cols;
        amgcl::backend::numa_vector<double> vals;
        long n_rows{0};
        long nnz{0};
    } Stencil;

public:

    Implementation(const CsrStencil &stencil, const std::vector<double> &mat_values,
                   double _relax_param, double _tol, unsigned _max_iter, bool _is_sym)
    : SystemStencil {
        .addr = std::vector<long>(stencil.addr().begin(), stencil.addr().end()),
        .cols = std::vector<long>(stencil.cols().begin(), stencil.cols().end()),
        .vals = mat_values,
        .n_rows = stencil.n_rows(),
        .nnz = stencil.n_nonzero()
    },
    RelaxParameter(_relax_param), Tolerance(_tol), MaxIterations(_max_iter), IsSystemSymmetric(_is_sym)
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

        if((*D_vals).size() != A.nrows)
            throw std::runtime_error("Solver init error: the main diagonal of system matrix contains zeroes, "
                                     "SOR(Symmetric SOR) method can not be used");

        (*D_vals).swap(DiagonalStencil.vals);

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

        GetTriangularComponent(UpperStencil);
        GetTriangularComponent(LowerStencil, false);

        InitMatrixWithoutOwning(U, UpperStencil);
        InitMatrixWithoutOwning(L, LowerStencil);

        /// Following parameter must stay in (0, 2) interval
        if(RelaxParameter <= 0.01 || RelaxParameter >= 1.99)
            SetOptionalRelaxParameter();

        if(!IsSystemSymmetric)
            CheckSystemForSymmetry();

    }

    [[nodiscard]] bool GetSymmetryStatus() const {
        return IsSystemSymmetric;
    }

    /**
     * @brief Following method implements solver procedure in Asymmetric system case.
     */
    void SolveSor() {};

    /**
     * @brief Following method implements solver procedure in Symmetric system case.
     */
    void SolveSSor() {};

private:

    /**
     * @brief Fills numa_vec with data from vec.
     * @tparam T
     * @param numa_vec
     * @param vec
     */
    template <typename T>
    void FillNumaFromVector(amgcl::backend::numa_vector<T>& numa_vec, const std::vector<T>& vec) {
        amgcl::backend::numa_vector<T> numa_tmp (vec);
        numa_vec.swap(numa_tmp);
    }

    /**
     * @brief This method implements checking for symmetry by scanning zeroes in [L + transpose(U)] matrix.
     */
    void CheckSystemForSymmetry()
    {
        auto transpose_upper_matrix_ptr = amgcl::backend::transpose(U);
        if(!transpose_upper_matrix_ptr)
            throw std::runtime_error("Check symmetry error: can not get transpose matrix");

        auto sum_result_matrix_ptr = amgcl::backend::sum(-1.0, L, 1.0, *transpose_upper_matrix_ptr);
        if(!sum_result_matrix_ptr)
            throw std::runtime_error("Check symmetry error: can not get sum of matrices");

        // TODO: nnz > real count of non zero elements...
        if(!sum_result_matrix_ptr->nnz)
            IsSystemSymmetric = true;
    }

    /**
    * @brief Method fills matrix with appropriate stencil data.
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
     * @brief Fills upper- / lower- triangular matrices stencil from incoming system stencil.
     * @param out_stencil
     * @param upper: upper triangular if true, lower triangular otherwise.
     */
    void GetTriangularComponent(Stencil& out_stencil, bool upper = true) {

        /// Creating upper triangular matrix U:
        std::vector<long> U_addr(SystemStencil.addr.size());
        std::vector<long> U_cols;
        std::vector<double> U_vals;

        // just fill with old (system matrices) values here
        for(auto i = 0; i < SystemStencil.addr.size(); ++i)
            U_addr[i] = SystemStencil.addr[i];

        for(auto i = 0; i < SystemStencil.addr.size() - 1; ++i)
        {
            for(auto j = SystemStencil.addr[i]; j < SystemStencil.addr[i + 1]; ++j)
            {
                if(upper)
                {
                    if(SystemStencil.cols[j] > i)
                    {
                        U_vals.emplace_back(SystemStencil.vals[j]);
                        U_cols.emplace_back(SystemStencil.cols[j]);
                    }
                    else
                        for(auto k = i; k < U_addr.size() - 1; ++k, U_addr[k]-=1);
                }
                else
                {
                    if(SystemStencil.cols[j] < i)
                    {
                        U_vals.emplace_back(SystemStencil.vals[j]);
                        U_cols.emplace_back(SystemStencil.cols[j]);
                    }
                    else
                        for(auto k = i; k < U_addr.size() - 1; ++k, U_addr[k]-=1);
                }
            }
        }

        FillNumaFromVector(out_stencil.addr, U_addr);
        FillNumaFromVector(out_stencil.cols, U_cols);
        FillNumaFromVector(out_stencil.vals, U_vals);

        out_stencil.n_rows = SystemStencil.n_rows;
        out_stencil.nnz = static_cast<long>(out_stencil.vals.size());
    }

    /**
     * @brief Implements the selection of optional relaxation parameter if input
     * value wasn't provided or was invalid.
     * @details w (aka omega) = 2 / 1 + sqrt(1 - spectral_radius*spectral_radius)\n
     * Or if w (aka omega) >= 2 after selection -> set w to UniversalRelaxParameter = 1.95
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

        /// Get upper bound of spectral radius by Gershgorin disk theorem with scaling by D^-1:
        double spectral_radius = amgcl::backend::spectral_radius<true, Matrix>(*product_matrix);

        RelaxParameter = 2.0 / (1.0 + sqrt(1.0 - spectral_radius * spectral_radius));

        if(RelaxParameter >= 2)
            RelaxParameter = UniversalRelaxParameter;

        /// Revert D matrices values
        D.val = DiagonalStencil.vals.data();
    }


//    /**
//     * @brief
//     */
//    void GetSchemasInverseComponent() {
//        if(!IsSystemSymmetric) {
//
//
//
//        }
//        else {
//
//
//
//        }
//    }


    Stencil SystemStencil, DiagonalStencil,
            LowerStencil, UpperStencil, SchemasInverseCompStencil;

    /// A = D + L + U factorization parts:
    Matrix D, L, U;

    /// Local solve parameters:
    /* Optional relaxation parameter: (0, 2); 1 - Gauss-Seidel case */
    double RelaxParameter;
    double Tolerance;
    unsigned MaxIterations;
    bool IsSystemSymmetric;

};

SorMatrixSolver::SorMatrixSolver(double RelaxParam, double Tol, unsigned MaxIter, bool IsSym)
: RelaxParameter(RelaxParam), Tolerance(Tol), MaxIterations(MaxIter), IsSystemSymmetric(IsSym)
{}

void SorMatrixSolver::set_matrix(const CsrStencil &Stencil, const std::vector<double> &MatValues)
{
    Impl = std::make_unique<SorMatrixSolver::Implementation>(Stencil, MatValues, RelaxParameter,
                                                             Tolerance, MaxIterations, IsSystemSymmetric);
    Impl->Initialize();
}

void SorMatrixSolver::solve(const std::vector<double> &rhs, std::vector<double> &ret) const {
    if(Impl)
    {
        if(Impl->GetSymmetryStatus())
            Impl->SolveSSor();
        else
        {
            throw std::runtime_error("TEST");
        }
    }
    else
        throw std::runtime_error("Solve error: matrix was not passed to the solver");
}