#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "common.hpp"
#include "appr/spatial_approximator.hpp"
#include "slae/a_matrix_solver.hpp"
#include "slae/sor_matrix_solver.hpp"

class PoissonSolver{
public:
	PoissonSolver(std::shared_ptr<ASpatialApproximator> appr);

	void set_bc_dirichlet(int btype, double value);
	void set_bc_dirichlet(int btype, std::function<double(Point)> value);

	void initialize();

	void solve(const std::function<double(Point)>& f, std::vector<double>& u);
	void solve(const std::vector<double>& f, std::vector<double>& u);
	void set_slae_solver(AMatrixSolver* solver);
private:
	std::shared_ptr<ASpatialApproximator> _approximator;

	std::vector<double> _mass;
	std::vector<double> _slae_rhs;

	std::unique_ptr<AMatrixSolver> _slae_solver;

	std::map<int, std::function<double(Point)>> _bc_dirichlet;
};

#endif
