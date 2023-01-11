#include "seidel_solver.hpp"
#include <iostream>
#include <math.h>

SeidelSolver::SeidelSolver(double eps, int max_iterations, int skip_res_iretations){
	this->eps = eps;
	this->max_iterations = max_iterations;
	this->skip_iterations = skip_res_iretations + 1;
}

void SeidelSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
	ret.resize(N);

	if (rhs.size() != N){
		throw std::runtime_error("Invalid Matrix-right side sizes");
	}

	std::vector<double> u(N, 0.0);

	int iter = 0;
	while (true){
		for (int i = 0; i < N; i++){
            double s = 0.0;
            for (int j = 1; j <= non_zeros[i]; j++) s -= v[i][j] * u[c[i][j]];
            s += rhs[i];
            u[i] = s / v[i][0];
        }
		
		if (iter % skip_iterations == 0 && solve_residual(u, rhs) < this->eps) break;
        
		iter += 1;
        if (iter >= this->max_iterations){
            std::cout << "WARN: The maximum number of iterations of " << this->max_iterations << " has been reached\n";
        	break;
        }
	}
	
	for(int i=0;i<N;i++) ret[i] = u[i];
}

void SeidelSolver::_check_matrix(int N) const{
    std::vector<int> bad_rows;
    for(int i=0; i < N ; i++){
        int ind1 = (this->stencil.addr())[i];
        int ind2 = (this->stencil.addr())[i+1]-1;

		double d = this->val[ind1];
        if (d <= 0) bad_rows.push_back(i);
		if (fabs(d) < 1e-15) throw std::runtime_error("Found zero on the diagonal");
    }

    if (bad_rows.size() > 0){
        std::cout << "WARN: SeidelSolver: Не выполнено достаточное условие сходимости в строках: ";
		for(int i=0; i < bad_rows.size(); i++){
			std::cout << bad_rows[i] << " ";
		}
		std::cout << "\n";
    }
}