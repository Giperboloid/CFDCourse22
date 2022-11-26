#ifndef MATRIX_SOLVER_HPP
#define MATRIX_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "a_matrix_solver.hpp"

// solve SLAE using amgc method
class AmgcMatrixSolver: public AMatrixSolver{
public:
	AmgcMatrixSolver(int maxit=1000, double eps=1e-8);
	~AmgcMatrixSolver();

	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;

	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const override;
private:
	int _maxit;
	double _tolerance;

	struct Impl;
	std::unique_ptr<Impl> _pimpl;
};

#endif
