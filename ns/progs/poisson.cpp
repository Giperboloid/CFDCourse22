#include <iostream>
#include "common.hpp"
#include "checks.hpp"
#include "prob/poisson_solver.hpp"
#include "grid/grid_builder.hpp"
#include "appr/linear_fem_approximator.hpp"

void linear_fem2(){
	// grid
	std::shared_ptr<Grid> grid = GridBuilder::build_from_gmshvtk("../../grids/rect1.vtk");

	// spatial approximator
	std::shared_ptr<LinearFemApproximator> linear_fem = LinearFemApproximator::build(grid);

	// solver
	PoissonSolver slv(linear_fem);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(4, 1);

	// rhs
	std::vector<double> rhs = linear_fem->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	linear_fem->vtk_save_scalar("u.vtk", x);

	// check solution
	CHECK_FLOAT3(x[150], 0);
	CHECK_FLOAT3(x[446], 0.477941);
	CHECK_FLOAT3(x[2], 1);
}

int main(){
	linear_fem2();
	std::cout << "DONE" << std::endl;
}
