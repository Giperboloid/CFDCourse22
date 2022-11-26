#ifndef TEST_SOLVER_HPP
#define TEST_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "a_matrix_solver.hpp"

// Тестовый класс. В качестве решения выдает линейную интерполяцию от 0.0 до значения value
class TestSolver: public AMatrixSolver{
public:
	TestSolver(double value=0.0);
	~TestSolver(){};

	void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) override;

	void solve(const std::vector<double>& rhs, std::vector<double>& ret) const override;
private:
	double value;
	CsrMatrix matrix;
};

#endif