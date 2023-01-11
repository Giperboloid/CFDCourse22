#include <iostream>
#include "a_jacobi_seidel_solver.hpp"

AJacobiSeidelSolver::AJacobiSeidelSolver() {}

void AJacobiSeidelSolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) {
	this->stencil = mat;
	this->val = mat_values;
	this->N = stencil.n_rows();

	this->c.resize(this->N);
	this->v.resize(this->N);
	this->non_zeros.resize(this->N);

	for (int i = 0; i < N; i++) {
		int ind1 = (this->stencil.addr())[i];
		int ind2 = (this->stencil.addr())[i + 1] - 1;
		c[i].resize(ind2 - ind1 + 1);
		v[i].resize(ind2 - ind1 + 1);
		for (int j = 0; j < c[i].size(); j++)
			c[i][j] = v[i][j] = 0.0;
		for (int j = ind1; j <= ind2; j++) {
			c[i][j - ind1] = (this->stencil.cols())[j];
			v[i][j - ind1] = (this->val)[j];
		}
		non_zeros[i] = ind2 - ind1;
	}

	this->_check_matrix(this->N);
}

double AJacobiSeidelSolver::solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const {
	std::vector<double> err(x.size(), 0.0);
	this->stencil.matvec(this->val, x, err);
	double s = 0.0;
	for (int i = 0;i < err.size(); i++) {
		s += (err[i] - rhs[i]) * (err[i] - rhs[i]);
	}
	return sqrt(s);
}