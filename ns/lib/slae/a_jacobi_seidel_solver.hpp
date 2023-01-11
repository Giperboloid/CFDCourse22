#ifndef A_JACOBI_SEIDEL_SOLVER_HPP
#define A_JACOBI_SEIDEL_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_matrix_solver.hpp"

class AJacobiSeidelSolver: public AMatrixSolver{
public:
	AJacobiSeidelSolver(double eps=0.001, int max_iterations=1000, int skip_res_iretations=0);

	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;
private:
	CsrStencil stencil;
	std::vector<double> val;
	double eps;
	int max_iterations;
	int skip_iterations;

	double solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const;
	virtual void _check_matrix(int N) = 0;

    std::vector<std::vector<int>> c;
	std::vector<std::vector<double>> v;
	std::vector<int> non_zeros;
	int N;
};

#endif