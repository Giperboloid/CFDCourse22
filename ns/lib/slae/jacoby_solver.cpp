#include "jacoby_solver.hpp"
#include <iostream>

JacobySolver::JacobySolver(double eps, int max_iterations, int skip_res_iretations){
	this->eps = eps;
	this->max_iterations = max_iterations;
	this->skip_iterations = skip_res_iretations + 1;
}

void JacobySolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values){
	this->stencil = mat;
	this->val = mat_values;
}

void JacobySolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
	int N = rhs.size();
	ret.resize(N);

	_check_Jacobi(N);

	std::vector<int> add = this->stencil.addr();
	std::vector<int> col = this->stencil.cols();
	std::vector<double> val = this->val;
	
	std::vector<double> u_old(N, 0.0);
	std::vector<double> u(N, 0.0);
	
	// Подготовка кэша матрицы
	std::vector<std::vector<int>> c;
	std::vector<std::vector<double>> v;
	std::vector<int> non_zeros(N, 0.0);
	c.resize(N);
	v.resize(N);
	for (int i = 0; i < N; i++){
		int ind1 = add[i];
        int ind2 = add[i+1]-1;
		c[i].resize(ind2-ind1+1);
		v[i].resize(ind2-ind1+1);
		for (int j = 0; j < c[i].size(); j++)
			c[i][j] = v[i][j] = 0.0;
		for (int j = ind1; j <= ind2; j++){
			c[i][j-ind1] = col[j];
			v[i][j-ind1] = val[j];
		}
		non_zeros[i] = ind2 - ind1;
	}

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


double JacobySolver::solve_residual(const std::vector<double>& x, const std::vector<double>& rhs) const{
	std::vector<double> err(x.size(), 0.0); 
	this->stencil.matvec(this->val, x, err);
	double s = 0.0;
	for (int i=0;i<err.size(); i++){
		s += (err[i] - rhs[i]) * (err[i] - rhs[i]);
	}
	return sqrt(s);
}

void JacobySolver::_check_Jacobi(int N) const{
	std::vector<int> add = this->stencil.addr();
	std::vector<int> col = this->stencil.cols();
	std::vector<double> val = this->val;
	
    std::vector<int> bad_rows;
    for(int i=0; i < N ; i++){
        int ind1 = add[i];
        int ind2 = add[i+1]-1;

		std::vector<double> v(ind2-ind1+1, 0.0);
		for (int j = ind1; j <= ind2; j++){
			v[j-ind1] = val[j];
		}
		double s = 0.0;
        for (int j = 1; j <= ind2 - ind1; j++){
			int vv = v[j];
			s += abs(vv);
		}
        if (s > abs(v[0])) bad_rows.push_back(i);
    }

    if (bad_rows.size() > 0){
        std::cout << "JacobySolver: Невыполнено достаточное условие сходимости в строках: ";
		for(int i=0; i < bad_rows.size(); i++){
			std::cout << bad_rows[i] << " ";
		}
		std::cout << "\n";
    }
}