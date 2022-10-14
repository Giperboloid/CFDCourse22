#include "psi_omega/psi_omega_semi_imp_vc.hpp"

class PsiOmega_SemiImpVC::VelocitySolver{
public:
	VelocitySolver(std::shared_ptr<AFemApproximator> appr):
			_approximator(appr.get()){
	}

	void solve(const std::vector<double>& psi, std::vector<double>& vx, std::vector<double>& vy){
		std::vector<Point> grad_psi = _approximator->element_centered_gradient(psi);
		for (int iel=0; iel<_approximator->n_elements(); ++iel){
			vx[iel] = grad_psi[iel].y;
			vy[iel] = -grad_psi[iel].x;
		}
	}
private:
	const AFemApproximator* _approximator;
};

class PsiOmega_SemiImpVC::OmegaSolver{
public:
	OmegaSolver(std::shared_ptr<AFemApproximator> appr):
			_approximator(appr.get()){
		_rhs.resize(_approximator->n_bases());
		_rhs_mat.resize(_approximator->stencil().n_nonzero());
		_slae_solver.reset(new AmgcMatrixSolver());
	}

	void set_tau(double tau){
		if (tau <= 0) throw std::runtime_error("Invalid tau value");
		_tau = tau;
	}

	double tau() const{
		return _tau;
	}

	void initialize(const APsiOmegaSolver* slv){
		// assemble vertices where omega = -Laplas(psi)
		std::set<int> bc_vert_set;
		for (auto it: slv->boundary_condition_by_class<PsiOmegaBc_Wall>()){
			for (int ivert: _approximator->boundary_bases(it.first)){
				bc_vert_set.insert(ivert);
			}
		}
		for (auto it: slv->boundary_condition_by_class<PsiOmegaBc_Input>()){
			for (int ivert: _approximator->boundary_bases(it.first)){
				bc_vert_set.insert(ivert);
			}
		}
		_bc_vert = std::vector<int>(bc_vert_set.begin(), bc_vert_set.end());

		// lhs = mass + tau / re * stiff
		const std::vector<double>& mass = _approximator->mass();
		const std::vector<double>& stiff = _approximator->stiff();
		std::vector<double> lhs(mass.size());
		for (size_t i=0; i<mass.size(); ++i){
			lhs[i] = mass[i] + _tau / slv->reynolds() * stiff[i];
		}
		// leave only mass for wall boundary
		const std::vector<int>& saddr = _approximator->stencil().addr();
		for (int ivert: _bc_vert){
			for (int a = saddr[ivert]; a<saddr[ivert+1]; ++a){
				lhs[a] = mass[a];
			}
		}

		// init solver
		_slae_solver->set_matrix(_approximator->stencil(), lhs);
	}

	void solve(const std::vector<double>& psi,
			const std::vector<double>& vx,
			const std::vector<double>& vy,
			const std::vector<double>& omega_old,
			std::vector<double>& omega){

		const CsrStencil& stencil = _approximator->stencil();
		const std::vector<double>& mass = _approximator->mass();
		const std::vector<double>& stiff = _approximator->stiff();

		// 1 assemble transport operator K(psi)
		std::vector<double> transport = _approximator->transport_v_elem(vx, vy, {});

		// 2 rhs = [Mass + tau * K(psi)] * w_old
		for (int i=0; i<stencil.n_nonzero(); ++i){
			_rhs_mat[i] = mass[i] + _tau * transport[i];
		}
		stencil.matvec(_rhs_mat, omega_old, _rhs);

		// 3 wall boundaries: rhs = stiff*psi
		for (int ivert: _bc_vert){
			_rhs[ivert] = stencil.matvec_irow(ivert, stiff, psi);
		}
		// 3.4 solve slae
		_slae_solver->solve(_rhs, omega);
	}
private:
	const AFemApproximator* _approximator;
	std::unique_ptr<AmgcMatrixSolver> _slae_solver;
	double _tau;

	std::vector<double> _rhs;
	std::vector<double> _rhs_mat;
	std::vector<int> _bc_vert;
};

PsiOmega_SemiImpVC::PsiOmega_SemiImpVC(std::shared_ptr<AFemApproximator> appr):
		APsiOmegaSolver(appr, appr),
		_approximator(appr.get()){
	
	_psi_solver.reset(new PoissonSolver(appr));
	_vel_solver.reset(new VelocitySolver(appr));
	_omega_solver.reset(new OmegaSolver(appr));
	_vx.resize(appr->n_elements());
	_vy.resize(appr->n_elements());
}

PsiOmega_SemiImpVC::~PsiOmega_SemiImpVC() = default;

void PsiOmega_SemiImpVC::set_tau(double tau){
	_omega_solver->set_tau(tau);
}

double PsiOmega_SemiImpVC::_compute_tau(){
	return _omega_solver->tau();
}

void PsiOmega_SemiImpVC::_solve_next_step(double tau){
	std::vector<double>& psi_old = _psi;
	std::vector<double>& omega_old = _omega;

	std::vector<double> psi(psi_old);
	std::vector<double> omega(omega_old);

	// 1. psi solution
	_psi_solver->solve(omega_old, psi);
	
	// 2. velocity solution
	_vel_solver->solve(psi, _vx, _vy);

	// 3. omega solution
	_omega_solver->solve(psi, _vx, _vy, omega_old, omega);

	// to the next layer
	std::swap(psi, psi_old);
	std::swap(omega, omega_old);
}

void PsiOmega_SemiImpVC::_initialize(){
	// 1. psi solver initialization
	// 1.1 boundary condition
	// input boundary conditions
	for (auto it: boundary_condition_by_class<PsiOmegaBc_Input>()){
		_psi_solver->set_bc_dirichlet(it.first, it.second->psifun);
	}
	// wall boundary conditions
	for (auto it: boundary_condition_by_class<PsiOmegaBc_Wall>()){
		_psi_solver->set_bc_dirichlet(it.first, it.second->psival);
	}
	// outflow boundary conditions:
	//    do nothing

	_psi_solver->initialize();
	_omega_solver->initialize(this);
}
