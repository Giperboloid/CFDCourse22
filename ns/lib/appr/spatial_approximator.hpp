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
	virtual void apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff) const = 0;
	virtual void apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> func, std::vector<double>& vec) const = 0;

	// =============== vtk data savers
	void vtk_save_scalar(std::string fname, const std::vector<double>& x, std::string dataname="data") const;
	void vtk_save_scalar(std::string fname, std::vector<const double*> x, std::vector<std::string> datanames) const;

	const Grid& grid() const { return *_grid; }
protected:
	virtual std::vector<double> _build_stiff() const = 0;
	virtual std::vector<double> _build_mass() const = 0;
	virtual CsrStencil _build_stencil() const = 0;

	ASpatialApproximator(std::shared_ptr<Grid> grid);

	std::shared_ptr<Grid> _grid;
private:
	struct Cache{
		std::vector<double> stiff;
		std::vector<double> mass;
		CsrStencil stencil;
		std::string vtk_outgrid;
	};
	mutable Cache _cache;

	// vtk output details
	virtual std::string _vtk_outgrid() const;
	virtual std::vector<double> _prepare_vtk_scalar(const double* x) const;
	virtual bool _vtk_use_cell_data() const;
};

#endif
