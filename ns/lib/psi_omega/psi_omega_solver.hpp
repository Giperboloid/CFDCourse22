#ifndef PSI_OMEGA_SOLVER_HPP
#define PSI_OMEGA_SOLVER_HPP

#include "common.hpp"
#include "prob/nonstationary_problem.hpp"
#include "appr/spatial_approximator.hpp"
#include "psi_omega/psi_omega_bc.hpp"

class APsiOmegaSolver: public ANonstationaryProblem{
public:
	APsiOmegaSolver(std::shared_ptr<ASpatialApproximator> psi_approx,
	                std::shared_ptr<ASpatialApproximator> omega_approx);
	virtual ~APsiOmegaSolver() = default;

	// setters
	void set_re(double re);
	void set_initial_condition_zero();
	template<typename T>
	void set_boundary_condition(int btype, T bc);

	// getters
	double re() const { return _re; }
	template<typename T>
	std::map<int, const T*> boundary_condition_by_class() const;

	const std::vector<double>& get_psi() const;
	const std::vector<double>& get_omega() const;
protected:
	double _re = 0;
	std::vector<double> _psi, _omega;
	std::shared_ptr<ASpatialApproximator> _psi_approximator;
	std::shared_ptr<ASpatialApproximator> _omega_approximator;
	std::map<int, std::unique_ptr<APsiOmegaBc>> _bc;

};

template<typename T>
void APsiOmegaSolver::set_boundary_condition(int btype, T bc){
	static_assert(std::is_base_of_v<APsiOmegaBc, T>);
	_bc[btype] = std::make_unique<T>(bc);
}

template<typename T>
std::map<int, const T*> APsiOmegaSolver::boundary_condition_by_class() const{
	std::map<int, const T*> ret;
	for (const auto& it: _bc){
		const T* b = dynamic_cast<const T*>(it.second.get());
		if (b != 0) ret.emplace(it.first, b);
	}
	return ret;
}

#endif
