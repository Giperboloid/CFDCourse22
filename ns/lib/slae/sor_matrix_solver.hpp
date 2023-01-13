#ifndef CFDAPP_SOR_MATRIX_SOLVER_HPP
#define CFDAPP_SOR_MATRIX_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_matrix_solver.hpp"

class SorMatrixSolver : public AMatrixSolver {

public:

    /**
     * @brief SOR solver constructor
     * @param RelaxParam - relaxation parameter (aka omega)
     * @param Tol - tolerance
     * @param MaxIter - upper bound of possible iterations count
     * @param IsSym - is SLAE symmetric
     */
    explicit SorMatrixSolver(double RelaxParam = 0.0, double Tol = 1e-8, unsigned MaxIter = 1000, bool IsSym = false);

    ~SorMatrixSolver() override = default;

    void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;
    void solve(const std::vector<double>& rhs, std::vector<double>& ret) const override;

private:

    class Implementation;
    std::unique_ptr<Implementation> Impl;

    double RelaxParameter;
    double Tolerance;
    unsigned MaxIterations;
    bool IsSystemSymmetric;
};



#endif //CFDAPP_SOR_MATRIX_SOLVER_HPP