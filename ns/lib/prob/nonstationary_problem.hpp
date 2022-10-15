#ifndef NONSTATIONARY_PROBLEM
#define NONSTATIONARY_PROBLEM

#include "common.hpp"
#include "prob/monitor.hpp"

// Abstract non stationary problem
class ANonstationaryProblem{
public:
	virtual ~ANonstationaryProblem() = default;

	// solves from tstart to tend
	void solve(double tend, double tstart=0);

	// adds solution monitor
	void add_monitor(std::shared_ptr<IMonitor> m);
protected:
	virtual double _compute_tau() = 0;
	virtual void _solve_next_step(double tau) = 0;
	virtual void _initialize();

	bool _initialized=false;
	std::vector<std::shared_ptr<IMonitor>> _monitors;
};


#endif
