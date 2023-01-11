#include "jacobi_solver.hpp"
#include <iostream>
#include <cmath>

JacobiSolver::JacobiSolver(double eps, int max_iterations, int skip_res_iretations) : AJacobiSeidelSolver(eps, max_iterations, skip_iterations){}

void JacobiSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const {
	if (rhs.size() != N) {
		throw std::runtime_error("Invalid Matrix-right side sizes");
	}

	if (ret.size() != N) {
		ret.resize(N);
		for (int i = 0; i < N; i++) ret[i] = 0.0;
	}

	std::vector<double> u_old(N, 0.0);

	bool convergence = false;
	for (int iter = 0; iter < this->max_iterations; iter++) {
		for (int i = 0; i < N; i++) {
			double s = 0.0;
			for (int j = 1; j <= non_zeros[i]; j++)
				s -= v[i][j] * u_old[c[i][j]];
			s += rhs[i];
			ret[i] = s / v[i][0];
		}

		if (iter % skip_iterations == 0 && solve_residual(ret, rhs) < this->eps) {
			convergence = true;
			break;
		}

		std::swap(u_old, ret);
	}
	if (!convergence) std::cout << "WARN: The maximum number of iterations of " << this->max_iterations << " has been reached\n";
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