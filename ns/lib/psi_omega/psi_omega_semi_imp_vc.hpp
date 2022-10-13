#ifndef PSI_OMEGA_SEMI_IMPLICIT_VC_HPP
#define PSI_OMEGA_SEMI_IMPLICIT_VC_HPP

#include "common.hpp"
#include "appr/fem_approximator.hpp"
#include "prob/poisson_solver.hpp"
#include "psi_omega/psi_omega_solver.hpp"

class PsiOmega_SemiImpVC: public APsiOmegaSolver{
public:
	PsiOmega_SemiImpVC(std::shared_ptr<AFemApproximator> appr);
	~PsiOmega_SemiImpVC();

	void set_tau(double tau);
protected:
	double _compute_tau() override;
	void _solve_next_step(double tau) override;
private:
	const AFemApproximator* _approximator;
	std::vector<double> _vx, _vy;

	std::unique_ptr<PoissonSolver> _psi_solver;

	class VelocitySolver;
	std::unique_ptr<VelocitySolver> _vel_solver;

	class OmegaSolver;
	std::unique_ptr<OmegaSolver> _omega_solver;

	void _initialize() override;
};

#endif
