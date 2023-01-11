#include "jacobi_solver.hpp"
#include <iostream>

JacobiSolver::JacobiSolver(double eps, int max_iterations, int skip_res_iretations)
{
	this->eps = eps;
	this->max_iterations = max_iterations;
	this->skip_iterations = skip_res_iretations + 1;
}

void JacobiSolver::solve(const std::vector<double> &rhs, std::vector<double> &ret) const
{
	int N = this->N;
	ret.resize(N);
	
	std::vector<double> u_old(N, 0.0);
	std::vector<double> u(N, 0.0);
	
	// Цикл решения
	int iter = 0;
	while (true){
		for (int i = 0; i < N; i++){          
			double s = 0.0;
			for (int j = 1; j <= non_zeros[i]; j++)
				s -= v[i][j] * u_old[c[i][j]];
			s += rhs[i];
			u[i] = s / v[i][0];
        }
		
		if (iter % skip_iterations == 0 && solve_residual(u, rhs) < this->eps) break;

		for (int i = 0; i < N; i++) u_old[i] = u[i];

		iter += 1;
		if (iter >= this->max_iterations){
			std::cout << "WARN: The maximum number of iterations of " << this->max_iterations << " has been reached\n";
			break;
		}
	}
	
	for(int i = 0; i < N; i++) ret[i] = u[i];
}

void JacobiSolver::_check_matrix(int N) const{
    std::vector<int> bad_rows;
    for(int i=0; i < N ; i++){
        int ind1 = (this->stencil.addr())[i];
        int ind2 = (this->stencil.addr())[i+1]-1;

		std::vector<double> v(ind2-ind1+1, 0.0);
		for (int j = ind1; j <= ind2; j++){
			v[j-ind1] = this->val[j];
		}
		double s = 0.0;
		for (int j = 1; j <= ind2 - ind1; j++){
			int vv = v[j];
			s += abs(vv);
		}
		if (s > abs(v[0])) bad_rows.push_back(i);
    }

	if (bad_rows.size() > 0){
		std::cout << "JacobiSolver: Невыполнено достаточное условие сходимости в строках: ";
		for(int i=0; i < bad_rows.size(); i++) std::cout << bad_rows[i] << " ";
		std::cout << "\n";
    }
}