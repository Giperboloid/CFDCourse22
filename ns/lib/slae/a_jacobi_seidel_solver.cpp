#include <iostream>
#include "a_jacobi_seidel_solver.hpp"

AJacobiSeidelSolver::AJacobiSeidelSolver(double eps, int max_iterations, int skip_res_iretations) {
	this->eps = eps;
	this->max_iterations = max_iterations;
	this->skip_iterations = skip_iterations + 2;
}

void AJacobiSeidelSolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) {
	this->stencil = mat;
	this->val = mat_values;
	this->N = mat.n_rows();

	this->c.resize(this->N);
	this->v.resize(this->N);
	this->non_zeros.resize(this->N);

	for (int i = 0; i < N; i++) {
		int ind1 = (mat.addr())[i];
		int ind2 = (mat.addr())[i + 1] - 1;
		c[i].resize(ind2 - ind1 + 1);
		v[i].resize(ind2 - ind1 + 1);
		for (int j = 0; j < c[i].size(); j++)
			c[i][j] = v[i][j] = 0.0;
		for (int j = ind1; j <= ind2; j++) {
			c[i][j - ind1] = (mat.cols())[j];
			v[i][j - ind1] = (this->val)[j];
		}
		non_zeros[i] = ind2 - ind1;
	}

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