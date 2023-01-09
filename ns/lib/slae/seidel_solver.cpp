#include "seidel_solver.hpp"
#include <iostream>

SeidelSolver::SeidelSolver(double eps, int max_iterations){
	this->eps = eps;
	this->max_iterations = max_iterations;
}

void SeidelSolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values){
	this->stencil = mat;
	this->val = mat_values;
}

void SeidelSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
	int N = rhs.size();
	ret.resize(N);

	_check_Seidel(N);

	std::vector<int> add = this->stencil.addr();
	std::vector<int> col = this->stencil.cols();
	std::vector<double> val = this->val;
	
	std::vector<double> u(N, 0.0);

	int iter = 0;
	while (true){
		for (int i = 0; i < N; i++){
            int ind1 = add[i];
            int ind2 = add[i+1]-1;

            std::vector<int> c(ind2-ind1+1, 0.0);
            std::vector<double> v(ind2-ind1+1, 0.0);
       		for (int j = ind1; j <= ind2; j++){
				c[j-ind1] = col[j];
				v[j-ind1] = val[j];
			}
            double s = 0.0;
            for (int j = 1; j <= ind2 - ind1; j++){
				int cc = c[j];
				int vv = v[j];
                s -= vv * u[cc];
            }
            s += rhs[i];
            u[i] = s / v[0];
        }
		
		if (solve_residual(u, rhs) < this->eps){
            break;
        }


		iter += 1;
        if (iter >= this->max_iterations){
            std::cout << "WARN: The maximum number of iterations of " << this->max_iterations << " has been reached\n";
        	break;
        }
	}
	
	for(int i=0;i<N;i++){
		ret[i] = u[i];
	}
}


double SeidelSolver::solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const{
	std::vector<double> err(x.size(), 0.0); 
	this->stencil.matvec(this->val, x, err);
	double s = 0.0;
	for (int i=0;i<err.size(); i++){
		s += (err[i] - rhs[i]) * (err[i] - rhs[i]);
	}
	return sqrt(s);
}

void SeidelSolver::_check_Seidel(int N) const{
	std::vector<int> add = this->stencil.addr();
	std::vector<int> col = this->stencil.cols();
	std::vector<double> val = this->val;
	
    std::vector<int> bad_rows;
    for(int i=0; i < N ; i++){
        int ind1 = add[i];
        int ind2 = add[i+1]-1;

		double d = val[ind1];
        if (d <= 0) bad_rows.push_back(i);
    }

    if (bad_rows.size() > 0){
        std::cout << "SeidelSolver: Невыполнено достаточное условие сходимости в строках: ";
		for(int i=0; i < bad_rows.size(); i++){
			std::cout << bad_rows[i] << " ";
		}
		std::cout << "\n";
    }
}