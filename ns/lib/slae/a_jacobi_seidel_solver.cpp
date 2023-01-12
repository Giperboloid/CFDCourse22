#include <iostream>
#include "a_jacobi_seidel_solver.hpp"

AJacobiSeidelSolver::AJacobiSeidelSolver(double eps, int max_iterations, int skip_res_iretations) {
	this->eps = eps;
	this->max_iterations = max_iterations;
	this->skip_iterations = skip_iterations + 1;
}

void AJacobiSeidelSolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) {
	this->stencil = mat;
	this->val = mat_values;
	this->N = mat.n_rows();

	this->non_zeros.resize(this->N);
	for (int i = 0; i < N; i++)
		non_zeros[i] = (mat.addr())[i + 1] - 1 - (mat.addr())[i];

	this->_check_matrix(this->N);
}

void AJacobiSeidelSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const
{
	if (rhs.size() != N) {
		throw std::runtime_error("Invalid Matrix-right side sizes");
	}

	if (ret.size() != N) {
		ret.resize(N);
		for (int i = 0; i < N; i++) ret[i] = 0.0;
	}

	if (!this->make_iterations(rhs, ret)) std::cout << "WARN: The maximum number of iterations of " << this->max_iterations << " has been reached\n";
}

double AJacobiSeidelSolver::solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const {
	std::vector<double> err(x.size(), 0.0);
	this->stencil.matvec(this->val, x, err);
	double s = 0.0;
	for (int i = 0; i < err.size(); i++) {
		s += (err[i] - rhs[i]) * (err[i] - rhs[i]);
	}
	return sqrt(s);
}