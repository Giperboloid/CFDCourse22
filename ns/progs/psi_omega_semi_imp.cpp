#include <iostream>
#include "common.hpp"
#include "geom.hpp"
#include "grid/grid_builder.hpp"
#include "psi_omega/psi_omega_semi_imp_vc.hpp"
#include "appr/linear_fem_approximator.hpp"

int main(){
	// === grid
	std::shared_ptr<Grid> grid = GridBuilder::build_from_gmshvtk(from_input_path("rect1.vtk"));

	// === spatial approximator
	std::shared_ptr<LinearFemApproximator> linear_fem = LinearFemApproximator::build(grid);

	// === solver
	PsiOmega_SemiImpVC slv(linear_fem);

	// === physical parameters
	slv.set_reynolds(100);

	// === solver parameters
	slv.set_tau(0.01);

	// ==== boundary conditions
	slv.set_boundary_condition(1, PsiOmegaBc_Input([](Point p){return p.y;}));
	slv.set_boundary_condition(2, PsiOmegaBc_Wall(0));
	slv.set_boundary_condition(3, PsiOmegaBc_Wall(1));

	// ==== initial conditions
	slv.set_initial_condition_zero();

	// ==== monitors
	// reports to console each 20 iterations
	auto console_monitor = std::make_shared<ConsoleIterReport>(20);
	slv.add_monitor(console_monitor);
	// saves psi and omega to vtk with dt = 1.0
	auto data_saver = std::make_shared<VtkFieldTimeSaver>(1, from_output_path("psi_omega"), linear_fem.get());
	data_saver->add_fun("psi", std::bind(&APsiOmegaSolver::psi, &slv));
	data_saver->add_fun("omega", std::bind(&APsiOmegaSolver::omega, &slv));
	slv.add_monitor(data_saver);

	// === solve for t=10
	slv.solve(10);

	std::cout << "DONE" << std::endl;
}
