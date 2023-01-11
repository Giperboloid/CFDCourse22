#ifndef A_JACOBI_SEIDEL_SOLVER_HPP
#define A_JACOBI_SEIDEL_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "slae/a_matrix_solver.hpp"

class AJacobiSeidelSolver : public AMatrixSolver {
public:
	AJacobiSeidelSolver(double eps, int max_iterations, int skip_res_iretations);
	virtual ~AJacobiSeidelSolver() = default;

	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;

protected:
	double solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const;

	virtual void _check_matrix(int N) const = 0;

	//Параметры решателя
	double eps;
	int max_iterations;
	int skip_iterations;

	// Описание матрицы
	CsrStencil stencil;
	std::vector<double> val;
	int N;

	//Кэш для ускорения работы с матрицей
	std::vector<std::vector<int>> c;
	std::vector<std::vector<double>> v;
	std::vector<int> non_zeros;
};

#endif