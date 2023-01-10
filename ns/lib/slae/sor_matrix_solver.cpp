#include "sor_matrix_solver.hpp"

struct SorMatrixSolver::Implementation {

public:

    Implementation(){};
    void solve();

};

SorMatrixSolver::SorMatrixSolver(unsigned MaxIter, double Tol)
: MaxIterations(MaxIter), Tolerance(Tol)
{}

void SorMatrixSolver::set_matrix(const CsrStencil &mat, const std::vector<double> &mat_values)
{}

void SorMatrixSolver::solve(const std::vector<double> &rhs, std::vector<double> &ret) const {
    if(Impl)
        Impl->solve();
    else
        throw std::runtime_error("Matrix was not passed to the solver");
}