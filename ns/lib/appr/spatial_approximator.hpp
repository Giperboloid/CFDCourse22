#ifndef SPATIAL_APPROXIMATOR_HPP
#define SPATIAL_APPROXIMATOR_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "grid/grid.hpp"

class ASpatialApproximator{
public:
	virtual ~ASpatialApproximator() = default;

	// =============== approximation dimension
	virtual int n_bases() const = 0;

	// =============== approximate function
	// approximates analytic function
	virtual std::vector<double> approximate(std::function<double(Point)> func) const = 0;

	// =============== stiff matrix stencil
	const CsrStencil& stencil() const;

	// =============== approximated differential operators
	const std::vector<double>& stiff() const;
	const std::vector<double>& mass() const;

	// =============== boundary conditions
	// get boundary basis indices with respective point coordinates
	const std::vector<std::pair<int, Point>>& boundary_bases(int ibnd) const;

	// forced dirichlet conditions
	void apply_bc_dirichlet(
		int ibnd,
		std::function<double(Point)> val_func,
		std::vector<double>& lhs, std::vector<double>& rhs) const;
	void apply_bc_dirichlet_lhs(int ibnd, std::vector<double>& lhs) const;
	void apply_bc_dirichlet_rhs(int ibnd, std::function<double(Point)> val_func, std::vector<double>& rhs) const;

	// stiff matrix conditions
	// du/dn = -q
	virtual void apply_bc_neumann_to_stiff(
		int ibnd,
		std::function<double(Point)> q_func,
		std::vector<double>& stiff, std::vector<double>& rhs) const;
	// du/dn = -alpha*u + beta
	virtual void apply_bc_robin_to_stiff(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::function<double(Point)> beta_func,
		std::vector<double>& stiff, std::vector<double>& rhs) const;

	// =============== vtk data savers
	void vtk_save_scalar(std::string fname, const std::vector<double>& x, std::string dataname="data") const;
	void vtk_save_scalar(std::string fname, std::vector<const double*> x, std::vector<std::string> datanames) const;

	const Grid& grid() const { return *_grid; }
protected:
	ASpatialApproximator(std::shared_ptr<Grid> grid);
	std::shared_ptr<Grid> _grid;
private:
	struct Cache{
		std::vector<double> stiff;
		std::vector<double> mass;
		CsrStencil stencil;
		std::map<int, std::vector<std::pair<int, Point>>> boundary_bases;
		std::string vtk_outgrid;
	};
	mutable Cache _cache;

	virtual std::vector<double> _build_stiff() const = 0;
	virtual std::vector<double> _build_mass() const = 0;
	virtual CsrStencil _build_stencil() const = 0;
	virtual std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const = 0;

	// vtk output details
	virtual std::string _vtk_outgrid() const;
	virtual std::vector<double> _prepare_vtk_scalar(const double* x) const;
	virtual bool _vtk_use_cell_data() const;
};

#endif
