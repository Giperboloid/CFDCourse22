#ifndef SEIDEL_SOLVER_HPP
#define SEIDEL_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_jacobi_seidel_solver.hpp"

class SeidelSolver : public AJacobiSeidelSolver {
public:
	SeidelSolver(double eps = 0.001, int max_iterations = 1000, int skip_res_iretations = 0);
	~SeidelSolver() {};

	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const override;
private:
	//Параметры решателя
	double eps;
	int max_iterations;
	int skip_iterations;

	void _check_matrix(int N) const override;
};

#endif