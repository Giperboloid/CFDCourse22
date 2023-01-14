
#include "sor_matrix_solver.hpp"

#include <amgcl/backend/builtin.hpp>
#include <numeric>

class SorMatrixSolver::Implementation {


public:
    Implementation(const CsrStencil& Stencil, const std::vector<double>& MatValues, double RelaxParameter,
        double Tolerance, unsigned MaxIterations, bool IsSystemSymmetric)
        :Matrix(Stencil), MatValues(MatValues), RelaxParameter(RelaxParameter), Tolerance(Tolerance), MaxIterations(MaxIterations), IsSystemSymmetric(IsSystemSymmetric)
    {}

    double addr_val(const std::vector<double> &data, const std::vector<int> &indices, const std::vector<int> &indptr, int irow, int jcol) {
        double val{ 0 };
        for (int a = indptr[irow]; a < indptr[irow + 1]; ++a) {
            if (indices[a] == jcol) {
                val = data[a];
            }
        }
        return val;
    }
    /**
     * @brief Following method implements solver procedure in Asymmetric system case.
     */
    void SolveSor(const std::vector<double>& rhs, std::vector<double>& ret)
    {
        
        int n = Matrix.n_rows();
        std::vector<double> xOld(n);
        std::vector<double> xNew(n);
        for (size_t i = 0; i < n; i++)
        {
            xOld[i] = 1.0;
            xNew[i] = 2.0;
        }
        int k = 0;
        double w = 1.95;
        double e = 10.0;
        std::vector<double> xkp1(n);
        std::vector<double> xk(n);

        std::vector<double> x(n);
        std::vector<double> A(n);

        std::vector<double> R(n);
        std::vector<double> A_1(n);
        std::vector<double> A_2(n);
        std::vector<double> d(n);
        for (size_t i = 0; i < n; i++)
        {
            d[i] = addr_val(MatValues, Matrix.cols(), Matrix.addr(), i, i);
        }
        while (k < MaxIterations && e > Tolerance)
        {
            for (size_t i = 0; i < n; i++)
            {
                //for (size_t j = 0; j < n; j++)
                //{
                //    if (j < i)
                //    {
                //        xkp1[j] = xNew[j];
                //    }
                //    else
                //    {
                //        xkp1[j] = 0.0;
                //    }
                //    if (j > i)
                //    {
                //        xk[j] = xOld[j];
                //    }
                //    else
                //    {
                //        xk[j] = 0.0;
                //    }
                //}
                //Matrix.matvec(MatValues, xkp1, A_1);
                //Matrix.matvec(MatValues, xk, A_2);
                //xNew[i] = (1.0 - w) * xOld[i] + w * (rhs[i] - A_1[i] - A_2[i]) / d[i];


                for (size_t j = 0; j < n; j++)
                {
                    if (j < i)
                    {
                        x[j] = xNew[j];
                    }
                    else if (j > i)
                    {
                        x[j] = xOld[j];
                    }
                    else if (i == j)
                    {
                        x[j] = 0.0;
                    }
                }
                Matrix.matvec(MatValues, x, A);
                xNew[i] = (1.0 - w) * xOld[i] + w * (rhs[i] - A[i]) / d[i];
            }
            xOld = xNew;
            k++;
            Matrix.matvec(MatValues, xNew, R);
            double R_sum {0};
            //for (size_t j = 0; j < n; j++)
            //{
            //    R[j] = (R[j] - rhs[j]) * (R[j] - rhs[j]);
            //    R_sum += R[j];
            //}

            std::transform(R.begin(), R.end(), rhs.begin(), R.begin(), std::minus<>());
            std::transform(R.begin(), R.end(), R.begin(), R.begin(), std::multiplies<>());
            R_sum = std::accumulate(R.begin(), R.end(),
                decltype(R)::value_type(0));

            e = sqrt(R_sum) / n;
            std::cout << "iter = "<< k << '\t' << "e = " << e << std::endl;
        }
        std::cout << "Solve" << std::endl;
        std::cout << "iter = "<< k << std::endl;
        ret.swap(xNew);
    };

    /**
     * @brief Following method implements solver procedure in Symmetric system case.
     */
    void SolveSSor() {};

private:

    const CsrStencil& Matrix;
    const std::vector<double> MatValues;
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
}

void SorMatrixSolver::solve(const std::vector<double> &rhs, std::vector<double> &ret) const {
    Impl->SolveSor(rhs, ret);
}