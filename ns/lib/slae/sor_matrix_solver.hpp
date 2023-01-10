#ifndef CFDAPP_SOR_MATRIX_SOLVER_HPP
#define CFDAPP_SOR_MATRIX_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_matrix_solver.hpp"

class SorMatrixSolver : public AMatrixSolver {

public:

    SorMatrixSolver(uint MaxIter = 1000, double Tol = 1e-8);

    ~SorMatrixSolver() override = default;

    void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;
    void solve(const std::vector<double>& rhs, std::vector<double>& ret) const override;

private:

    struct Implementation;
    std::unique_ptr<Implementation> Impl;

    double Tolerance;
    uint MaxIterations;

    uint8_t RelaxParam {1};
};



#endif //CFDAPP_SOR_MATRIX_SOLVER_HPP