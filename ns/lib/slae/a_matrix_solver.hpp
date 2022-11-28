#ifndef A_MATRIX_SOLVER_HPP
#define A_MATRIX_SOLVER_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"

class AMatrixSolver{
public:
	AMatrixSolver(){};
	virtual ~AMatrixSolver() = default;

	void set_matrix(const CsrMatrix& mat);
	virtual void set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values) = 0;

	virtual void solve(const std::vector<double>& rhs, std::vector<double>& ret) const = 0;	
};

#endif
