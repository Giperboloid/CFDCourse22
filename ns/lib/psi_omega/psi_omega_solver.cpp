#include "psi_omega/psi_omega_solver.hpp"

APsiOmegaSolver::APsiOmegaSolver(std::shared_ptr<ASpatialApproximator> psi_approx,
                                 std::shared_ptr<ASpatialApproximator> omega_approx):
		_psi_approximator(psi_approx),
		_omega_approximator(omega_approx){
	_psi.resize(_psi_approximator->n_bases());
	_omega.resize(_omega_approximator->n_bases());
}

void APsiOmegaSolver::set_re(double re){
	if (re <= 0) throw std::runtime_error("Invalid Reynolds number");
	_re = re;
}

void APsiOmegaSolver::set_initial_condition_zero(){
	std::fill(_psi.begin(), _psi.end(), 0.0);
	std::fill(_omega.begin(), _omega.end(), 0.0);
}

const std::vector<double>& APsiOmegaSolver::psi() const{
	return _psi;
}
const std::vector<double>& APsiOmegaSolver::omega() const{
	return _omega;
}
