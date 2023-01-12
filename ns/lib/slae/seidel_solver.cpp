#include "seidel_solver.hpp"
#include <iostream>
#include <cmath>

SeidelSolver::SeidelSolver(double eps, int max_iterations, int skip_res_iretations): AJacobiSeidelSolver(eps, max_iterations, skip_res_iretations) {}

void SeidelSolver::make_iteration(const std::vector<double>& rhs, std::vector<double>& ret) const {
	for (int i = 0; i < N; i++) {
		int ind1 = (this->stencil.addr())[i];
		double s = 0.0;
		for (int j = 1; j <= non_zeros[i]; j++)
			s -= (this->val)[j + ind1] * ret[(this->stencil.cols())[j + ind1]];
		s += rhs[i];
		ret[i] = s / (this->val)[ind1];
	}
}

void SeidelSolver::create_solver_cache() const {}

void SeidelSolver::_check_matrix(int N) const {
	std::vector<int> bad_rows;
	for (int i = 0; i < N; i++) {
		int ind1 = (this->stencil.addr())[i];
		int ind2 = (this->stencil.addr())[i + 1] - 1;

		double d = this->val[ind1];
		if (d <= 0) bad_rows.push_back(i);
		if (std::fabs(d) < 1e-15) throw std::runtime_error("Found zero on the diagonal");
	}

	if (bad_rows.size() > 0) {
		std::cout << "WARN: SeidelSolver: Не выполнено достаточное условие сходимости в строках: ";
		for (int i = 0; i < bad_rows.size(); i++) {
			std::cout << bad_rows[i] << " ";
		}
		std::cout << "\n";
	}
}