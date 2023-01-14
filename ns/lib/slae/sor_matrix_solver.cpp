
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

        DiagonalStencil.nnz = static_cast<long>(DiagonalStencil.vals.size());
        DiagonalStencil.cols.resize(DiagonalStencil.nnz);

        DiagonalStencil.n_rows = static_cast<long>(A.nrows);
        DiagonalStencil.addr.resize(DiagonalStencil.n_rows + 1);


        for(long i = 0; i <= DiagonalStencil.n_rows; ++i) {
            if (i != DiagonalStencil.n_rows) {
                DiagonalStencil.addr[i] = i;
                DiagonalStencil.cols[i] = i;
            } else
                DiagonalStencil.addr[i] = i;
        }

        InitMatrixWithoutOwning(D, DiagonalStencil);

        GetTriangularComponent(UpperStencil);
        GetTriangularComponent(LowerStencil, false);

        InitMatrixWithoutOwning(U, UpperStencil);
        InitMatrixWithoutOwning(L, LowerStencil);

        /// Following parameter must stay in (0, 2) interval
        if(RelaxParameter <= 0.0 || RelaxParameter >= 2.0)
            SetOptionalRelaxParameter();

        if(!IsSystemSymmetric)
            CheckSystemForSymmetry();

        GetSchemasInverseComponent();

    }

    [[nodiscard]] bool GetSymmetryStatus() const {
        return IsSystemSymmetric;
    }

    /**
     * @brief Following method implements solver procedure in Asymmetric system case.
     * @param rhs
     * @param ret
     */
    void SolveSor(const std::vector<double>& rhs, std::vector<double>& ret) const {

        using amgcl::backend::sum, amgcl::backend::product;

        /// First approximation - unit vector {1, 1, 1, 1, ...}
        std::vector<double> X (SystemStencil.n_rows, 1);

        /// Convert X and Rhs vectors to matrices:
        Matrix XMatrix, RhsMatrix;

        Stencil XStencil{.vals = X, .n_rows = SystemStencil.n_rows, .nnz = SystemStencil.n_rows};
        XStencil.cols.resize(SystemStencil.n_rows);
        XStencil.addr.resize(SystemStencil.n_rows + 1);

        Stencil RhsStencil{.vals = rhs, .n_rows = SystemStencil.n_rows, .nnz = SystemStencil.n_rows};
        RhsStencil.cols.resize(SystemStencil.n_rows);
        RhsStencil.addr.resize(SystemStencil.n_rows + 1);

        for(auto i = 0; i < SystemStencil.n_rows + 1; ++i) {
                XStencil.addr[i] = i;
                RhsStencil.addr[i] = i;
        }

        InitMatrixWithoutOwning(XMatrix, XStencil);

        EraseZeroes(RhsStencil);
        InitMatrixWithoutOwning(RhsMatrix, RhsStencil);


        auto XMultiplierMatrix = sum(RelaxParameter, U, RelaxParameter - 1, D);
        if(!XMultiplierMatrix)
            throw std::runtime_error("Solve SOR error: can not get matrices sum [wU + (w - 1)D]");

        /// Solve
        size_t iter_count = 0;
        while(++iter_count <= MaxIterations)
        {
            auto XMultiplier_by_ = product(*XMultiplierMatrix, XMatrix);
            if(!XMultiplier_by_)
                throw std::runtime_error("Solve SOR error: can not get product [wU + (w - 1)D].X");

            auto WRhs_minus_ = sum(RelaxParameter, RhsMatrix, -1.0, *XMultiplier_by_);
            if(!WRhs_minus_)
                throw std::runtime_error("Solve SOR error: can not get sum of wB - [wU + (w - 1)D].X");

            auto InverseComp_by_ = product(InverseComponent, *WRhs_minus_);
            if(!InverseComp_by_)
                throw std::runtime_error("Solve SOR error: can not get product of [D + wL]^-1 . [wB - [wU + (w - 1)D].X]");

            XMatrix = *InverseComp_by_;

            for(auto i = 0; i < XMatrix.nnz; ++i)
                std::cout << XMatrix.val[i] << std::endl;

            std::cout << std::endl << std::endl;
        }

        for(auto i = 0; i < X.size(); ++i)
            X[i] = XMatrix.val[i];

        ret.swap(X);
    };

    /**
     * @brief Following method implements solver procedure in Symmetric system case.
     * @param Rhs
     * @param ret
     */
    void SolveSSor(const std::vector<double>& Rhs, std::vector<double>& ret) const {
        throw std::runtime_error("Solve error: S-SOR solver is not implemented");
    };

private:

    /**
     * @brief Fills numa_vec with data from vec.
     * @tparam T
     * @param numa_vec
     * @param vec
     */
    template <typename T>
    static void FillNumaValsFromVector(amgcl::backend::numa_vector<T>& numa_vec, const std::vector<T>& vec) {
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
            throw std::runtime_error("Check symmetry error: can not get transpose matrix U^T");

        auto sum_result_matrix_ptr = amgcl::backend::sum(-1.0, L, 1.0, *transpose_upper_matrix_ptr);
        if(!sum_result_matrix_ptr)
            throw std::runtime_error("Check symmetry error: can not get sum of matrices -L + U^T");

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
        std::vector<long> U_addr(SystemStencil.addr.size()), U_cols;
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

        FillNumaValsFromVector(out_stencil.addr, U_addr);
        FillNumaValsFromVector(out_stencil.cols, U_cols);
        FillNumaValsFromVector(out_stencil.vals, U_vals);

        out_stencil.n_rows = SystemStencil.n_rows;
        out_stencil.nnz = static_cast<long>(out_stencil.vals.size());

        EraseZeroes(out_stencil);
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
            throw std::runtime_error("Setting relaxation parameter error: can not get sum of triangle matrices U + L");

        amgcl::backend::numa_vector<double>Inverse_D_vals(DiagonalStencil.vals.size());

        for(auto i = 0; i < Inverse_D_vals.size(); ++i)
            Inverse_D_vals[i] = 1.0 / DiagonalStencil.vals[i];

        /// Temporary set D matrices values pointer to 'Inverse_D_vals' <-> gives opportunity to use D as D^-1
        /// without creating new one matrix
        D.val = Inverse_D_vals.data();

        auto product_matrix = amgcl::backend::product(D, *triangles_sum_ptr);
        if(!product_matrix)
            throw std::runtime_error("Setting relaxation parameter error: can not get matrices product (D^-1)(U + L)");

        /// Get upper bound of spectral radius by Gershgorin disk theorem with scaling by D^-1:
        double spectral_radius = amgcl::backend::spectral_radius<true, Matrix>(*product_matrix);

        RelaxParameter = 2.0 / (1.0 + sqrt(1.0 - spectral_radius * spectral_radius));

        if(RelaxParameter >= 2)
            RelaxParameter = UniversalRelaxParameter;

        /// Revert D matrices values
        D.val = DiagonalStencil.vals.data();
    }

    /**
     * @brief
     * @param stencil
     */
    static void EraseZeroes(Stencil& stencil) {

        std::vector<double> vals_new;
        std::vector<long> addr_new(stencil.addr.size()), cols_new;

        for (auto i = 0; i < stencil.addr.size(); ++i)
            addr_new[i] = stencil.addr[i];

        for (auto i = 0; i < stencil.addr.size() - 1; ++i)
        {
            for (auto j = stencil.addr[i]; j < stencil.addr[i + 1]; ++j)
            {
                if (stencil.vals[j] != 0.0)
                {
                    vals_new.emplace_back(stencil.vals[j]);
                    cols_new.emplace_back(stencil.cols[j]);
                }
                else
                    for (int k = i; k < addr_new.size() - 1; ++k, addr_new[k] -= 1);
            }
        }

        FillNumaValsFromVector(stencil.addr, addr_new);
        FillNumaValsFromVector(stencil.cols, cols_new);
        FillNumaValsFromVector(stencil.vals, vals_new);

        stencil.nnz = static_cast<long>(vals_new.size());
    }

    /**
     * @brief
     * @param val
     * @param col
     * @param ptr
     * @param row_idx
     * @param col_idx
     * @return
     */
    static double GetMatrixValue (const double* val, const long* col, const long* ptr,
                                            long row_idx, long col_idx) {
        double return_val{ 0 };
        for (auto a = ptr[row_idx]; a < ptr[row_idx + 1]; ++a) {
            if (col[a] == col_idx) {
                return_val = val[a];
                break;
            }
        }
        return return_val;
    }

    /**
     * @brief
     */
    void GetSchemasInverseComponent() {

        if(!IsSystemSymmetric) {

            /// Schema's inverse component is: (D + wL)^-1;
            auto inverse_comp_ptr = amgcl::backend::sum(1.0, D, RelaxParameter, L);
            if(!inverse_comp_ptr)
                throw std::runtime_error("Getting schema's inverse component error: can not get the sum of matrices [D + wL]");

            std::vector<double> vals_new;
            std::vector<long> cols_new(inverse_comp_ptr->nnz), addr_new(inverse_comp_ptr->nrows + 1);

            for(auto i = 0; i < cols_new.size(); ++i)
                cols_new[i] = inverse_comp_ptr->col[i];

            for(auto i = 0; i < addr_new.size(); ++i)
                addr_new[i] = inverse_comp_ptr->ptr[i];

            /// Fill inverse component's diagonal elements:
            for(auto i = 0; i < inverse_comp_ptr->nrows; ++i) {
                for(auto j = inverse_comp_ptr->ptr[i]; j < inverse_comp_ptr->ptr[i + 1]; ++j){
                    if(inverse_comp_ptr->col[j] == i)
                        vals_new.emplace_back(1.0 / inverse_comp_ptr->val[j]);
                    else
                        vals_new.emplace_back(inverse_comp_ptr->val[j]);
                }
            }

            /// Fill below diagonal elements:
            for(auto i = 0; i < inverse_comp_ptr->nrows; ++i) {
                for(auto j = inverse_comp_ptr->ptr[i]; j < inverse_comp_ptr->ptr[i + 1]; ++j) {
                    if(inverse_comp_ptr->col[j] < i)
                    {
                        double sum { 0.0 };
                        for (auto k = 0; k < i; ++k) {
                            sum += GetMatrixValue(vals_new.data(), cols_new.data(), addr_new.data(), k, cols_new[j])
                                   * GetMatrixValue(inverse_comp_ptr->val, cols_new.data(), addr_new.data(), i, k);
                        }

                        vals_new[j] = -sum / GetMatrixValue(inverse_comp_ptr->val, cols_new.data(), addr_new.data(), i, i);
                    }
                }
            }

            FillNumaValsFromVector(SchemasInverseCompStencil.addr, addr_new);
            FillNumaValsFromVector(SchemasInverseCompStencil.cols, cols_new);
            FillNumaValsFromVector(SchemasInverseCompStencil.vals, vals_new);

            SchemasInverseCompStencil.n_rows = static_cast<long>(inverse_comp_ptr->nrows);
            SchemasInverseCompStencil.nnz = static_cast<long>(vals_new.size());
        }

        EraseZeroes(SchemasInverseCompStencil);
        InitMatrixWithoutOwning(InverseComponent, SchemasInverseCompStencil);
    }


    Stencil SystemStencil, DiagonalStencil,
            LowerStencil, UpperStencil, SchemasInverseCompStencil;

    /// A = D + L + U factorization parts & schema's inverse component:
    Matrix D, L, U, InverseComponent;

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

void SorMatrixSolver::solve(const std::vector<double> &Rhs, std::vector<double> &ret) const {
    if(Impl)
    {
        if(Impl->GetSymmetryStatus())
            Impl->SolveSSor(Rhs, ret);
        else
            Impl->SolveSor(Rhs, ret);
    }
    else
        throw std::runtime_error("Solve error: matrix was not passed to the solver");
}