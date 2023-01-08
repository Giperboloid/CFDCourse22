#include "jacoby_solver.hpp"
#include <iostream>

JacobySolver::JacobySolver(double eps, int max_iterations){
	this->eps = eps;
	this->max_iterations = max_iterations;
}

void JacobySolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values){
	this->stencil = mat;
	this->val = mat_values;
}

void JacobySolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
	int N = rhs.size();
	ret.resize(N);

	std::vector<int> add = this->stencil.addr();
	std::vector<int> col = this->stencil.cols();
	std::vector<double> val = this->val;
	
	std::vector<double> u_old(N, 0.0);
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
                s -= vv * u_old[cc];
            }
            s += rhs[i];
            u[i] = s / v[0];
        }
		
		if (solve_residual(u, rhs) < this->eps){
            break;
        }

		for (int i=0; i< N; i++){
			u_old[i] = u[i];
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


double JacobySolver::solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const{
	std::vector<double> err(x.size(), 0.0); 
	this->stencil.matvec(this->val, x, err);
	double s = 0.0;
	for (int i=0;i<err.size(); i++){
		s += (err[i] - rhs[i]) * (err[i] - rhs[i]);
	}
	return sqrt(s);
}