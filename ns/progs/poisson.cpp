#include <iostream>
#include <string>
#include <time.h>

#include "common.hpp"
#include "prog_common.hpp"
#include "prob/poisson_solver.hpp"
#include "grid/grid_builder.hpp"
#include "appr/linear_fem_approximator.hpp"
#include "appr/fvm_approximator.hpp"

#include "slae/a_matrix_solver.hpp"
#include "slae/a_jacobi_seidel_solver.hpp"
#include "slae/test_solver.hpp"
#include "slae/jacobi_solver.hpp"
#include "slae/seidel_solver.hpp"
#include "slae/amgcl_matrix_solver.hpp"

void linear_fvm1(AMatrixSolver* solver){
	// grid
	std::shared_ptr<Grid> grid = GridBuilder::build_regular1(10, 101);

	// spatial approximator
	std::shared_ptr<FvmApproximator> appr = FvmApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);

	// пример установки другого решателя
	slv.set_slae_solver(solver);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(2, 1);

	// rhs
	std::vector<double> rhs = appr->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fvm.vtk"), x);
	
	// check solution
	CHECK_FLOAT3(x[0], 0.005);
	CHECK_FLOAT3(x[5], 0.055);
	CHECK_FLOAT3(x[99], 0.995);
}

void linear_fvm2(AMatrixSolver* solver){
	// grid
	std::string grid_filename = from_input_path("rect1.vtk");
	std::shared_ptr<Grid> grid = GridBuilder::build_from_gmshvtk(grid_filename);

	// spatial approximator
	std::shared_ptr<FvmApproximator> appr = FvmApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);
	slv.set_slae_solver(solver);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(4, 1);

	// rhs
	std::vector<double> rhs = appr->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fvm.vtk"), x);

	// check solution
	CHECK_FLOAT3(x[819], 0.148651);
}

void linear_fem2(AMatrixSolver* solver){
	// grid
	std::string grid_filename = from_input_path("rect1.vtk");
	std::shared_ptr<Grid> grid = GridBuilder::build_from_gmshvtk(grid_filename);

	// spatial approximator
	std::shared_ptr<LinearFemApproximator> appr = LinearFemApproximator::build(grid);

	// solver
	PoissonSolver slv(appr);
	slv.set_slae_solver(solver);

	// bc
	slv.set_bc_dirichlet(1, 0);
	slv.set_bc_dirichlet(4, 1);

	// rhs
	std::vector<double> rhs = appr->approximate([](Point p){ return 0; });

	// solve
	std::vector<double> x;
	slv.initialize();
	slv.solve(rhs, x);

	// show solution
	appr->vtk_save_scalar(from_output_path("poisson_fem.vtk"), x);

	// check solution
	CHECK_FLOAT3(x[150], 0);
	CHECK_FLOAT3(x[446], 0.477941);
	CHECK_FLOAT3(x[2], 1);
}

int main(){
	try{
		std::cout << "\n\nAmgcMatrixSolver\n\n";
		linear_fvm1(new AmgcMatrixSolver());
		linear_fvm2(new AmgcMatrixSolver());
		linear_fem2(new AmgcMatrixSolver());

		std::cout << "\n\nJacobySolver\n\n";
		linear_fvm1(new JacobiSolver(0.0001, 100000, 100));
		linear_fvm2(new JacobiSolver(0.00001, 100000, 10000));
		linear_fem2(new JacobiSolver(0.00001, 100000, 10000));

		std::cout << "\n\nSeidelSolver\n\n";
		linear_fvm1(new SeidelSolver(0.0001, 100000, 100));
		linear_fvm2(new SeidelSolver(0.00001, 100000, 10000));
		linear_fem2(new SeidelSolver(0.00001, 100000, 10000));

		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << " " << e.what() << std::endl;
	}
}
