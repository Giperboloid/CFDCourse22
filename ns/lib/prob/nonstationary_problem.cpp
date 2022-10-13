#include "prob/nonstationary_problem.hpp"

void ANonstationaryProblem::solve(double tend, double tstart){
	if (!_initialized){
		std::cout << "=== Initialization" << std::endl;
		_initialize();
		for (auto& m: _monitors){
			m->initialize(tstart);
		}
		_initialized=true;
	}

	std::cout << "=== Start computation at t = " << tstart << std::endl;

	double t = tstart;
	while (t < tend - 1e-8){
		double tau = _compute_tau();
		if (t + tau > tend){
			tau = tend - t;
		}
		_solve_next_step(tau);
		t += tau;

		bool forced_break = false;
		for (auto& m: _monitors){
			if (m->run(t)) forced_break = true;
		}

		if (forced_break) break;
	}

	for (auto& m: _monitors){
		m->finalize(t);
	}

	std::cout << "=== Stop computation at t = " << t << std::endl;
}

void ANonstationaryProblem::add_monitor(std::shared_ptr<IMonitor> m){
	_monitors.push_back(m);
}

void ANonstationaryProblem::_initialize(){
}
