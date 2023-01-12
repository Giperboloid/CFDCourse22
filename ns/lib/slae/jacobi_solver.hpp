#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_jacobi_seidel_solver.hpp"

class JacobiSolver : public AJacobiSeidelSolver {
public:
	JacobiSolver(double eps = 0.001, int max_iterations = 1000, int skip_res_iretations = 0);
	~JacobiSolver() {};

private:
	void make_iteration(const std::vector<double>& rhs, std::vector<double>& ret) const override;
	void create_solver_cache() const override;
	void _check_matrix(int N) const override;

	std::vector<double>* u_old = new std::vector<double>;
};

#endif