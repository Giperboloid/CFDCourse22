#include "seidel_solver.hpp"
#include <iostream>
#include <cmath>

SeidelSolver::SeidelSolver(double eps, int max_iterations, int skip_res_iretations): AJacobiSeidelSolver(eps, max_iterations, skip_iterations) {}

bool SeidelSolver::make_iterations(const std::vector<double>& rhs, std::vector<double>& ret) const {
	bool convergence = false;
	for (int iter = 0; iter < this->max_iterations; iter++) {
		for (int i = 0; i < N; i++) {
			double s = 0.0;
			for (int j = 1; j <= non_zeros[i]; j++) s -= v[i][j] * ret[c[i][j]];
			s += rhs[i];
			ret[i] = s / v[i][0];
		}

		if (iter % skip_iterations == 0 && solve_residual(ret, rhs) < this->eps) {
			convergence = true;
			break;
		}
	}
	return convergence;
}

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