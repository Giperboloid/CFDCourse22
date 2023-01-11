#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_matrix_solver.hpp"

class JacobiSolver: public AMatrixSolver{
public:
	JacobiSolver(double eps=0.001, int max_iterations=1000, int skip_res_iretations=0);
	~JacobiSolver(){};

	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;

	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const override;
private:
	CsrStencil stencil;
	std::vector<double> val;
	double eps;
	int max_iterations;
	int skip_iterations;

	double solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const;
	void _check_Jacobi(int N) const;
};

#endif