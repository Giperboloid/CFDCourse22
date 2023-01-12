#include "jacobi_solver.hpp"
#include <iostream>
#include <cmath>

JacobiSolver::JacobiSolver(double eps, int max_iterations, int skip_res_iretations): AJacobiSeidelSolver(eps, max_iterations, skip_res_iretations) {}

void JacobiSolver::make_iteration(const std::vector<double>& rhs, std::vector<double>& ret) const {
	std::swap(this->u_old, ret);
	for (int i = 0; i < N; i++) {
		int ind1 = (this->stencil.addr())[i];
		double s = 0.0;
		for (int j = 1; j <= non_zeros[i]; j++)
			s -= (this->val)[j + ind1] * this->u_old[(this->stencil.cols())[j + ind1]];

		s += rhs[i];
		ret[i] = s / (this->val)[ind1];
	}
}

void JacobiSolver::create_solver_cache() const {
	this->u_old.resize(this->N);
}

void JacobiSolver::_check_matrix(int N) const {
	std::vector<int> bad_rows;
	for (int i = 0; i < N; i++) {
		int ind1 = (this->stencil.addr())[i];
		int ind2 = (this->stencil.addr())[i + 1] - 1;

		std::vector<double> v(ind2 - ind1 + 1, 0.0);
		for (int j = ind1; j <= ind2; j++) {
			v[j - ind1] = this->val[j];
		}
		double s = 0.0;
		for (int j = 1; j <= ind2 - ind1; j++) {
			int vv = v[j];
			s += std::fabs(vv);
		}
		if (s > std::fabs(v[0])) bad_rows.push_back(i);
		if (std::fabs(v[0]) < 1e-15) throw std::runtime_error("Found zero on the diagonal");
	}

	if (bad_rows.size() > 0) {
		std::cout << "WARN: JacobiSolver: Не выполнено достаточное условие сходимости в строках: ";
		for (int i = 0; i < bad_rows.size(); i++) std::cout << bad_rows[i] << " ";
		std::cout << "\n";
	}
}