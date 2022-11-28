#include "a_matrix_solver.hpp"

void AMatrixSolver::set_matrix(const CsrMatrix& mat){
	set_matrix(mat, mat.vals());
}
