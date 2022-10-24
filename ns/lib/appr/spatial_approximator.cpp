#include "appr/spatial_approximator.hpp"
#include <fstream>


ASpatialApproximator::ASpatialApproximator(std::shared_ptr<Grid> grid){
	_grid = grid;
}

const CsrStencil& ASpatialApproximator::stencil() const{
	if (_cache.stencil.n_rows() == 0){
		_cache.stencil = _build_stencil();
	}
	return _cache.stencil;
}
const std::vector<double>& ASpatialApproximator::stiff() const{
	if (_cache.stiff.empty()){
		_cache.stiff = _build_stiff();
	}
	return _cache.stiff;
}
const std::vector<double>& ASpatialApproximator::mass() const{
	if (_cache.mass.empty()){
		_cache.mass = _build_mass();
	}
	return _cache.mass;
}

void ASpatialApproximator::apply_bc_dirichlet(
		int ibnd,
		std::function<double(Point)> val_func,
		std::vector<double>& lhs, std::vector<double>& rhs) const{
	const CsrStencil& s = stencil();

	const std::vector<std::pair<int, Point>>& bnd_points = boundary_bases(ibnd);

	if (!lhs.empty())
	for (std::pair<int, Point> it: bnd_points){
		s.set_unit_diagonal(it.first, lhs);
	}

	if (!rhs.empty())
	for (std::pair<int, Point> it: bnd_points){
		rhs[it.first] = val_func(it.second);
	}
}

void ASpatialApproximator::apply_bc_dirichlet_lhs(int ibnd, std::vector<double>& lhs) const{
	std::vector<double> rhs;
	return apply_bc_dirichlet(ibnd, [](Point){ return 0; }, lhs, rhs);
}

void ASpatialApproximator::apply_bc_dirichlet_rhs(int ibnd, std::function<double(Point)> val_func, std::vector<double>& rhs) const{
	std::vector<double> lhs;
	return apply_bc_dirichlet(ibnd, val_func, lhs, rhs);
}

const std::vector<std::pair<int, Point>>& ASpatialApproximator::boundary_bases(int ibnd) const{
	if (_cache.boundary_bases.empty()){
		_cache.boundary_bases = _build_boundary_bases();
	}
	return _cache.boundary_bases[ibnd];
}

void ASpatialApproximator::apply_bc_neumann_to_stiff(
		int ibnd,
		std::function<double(Point)> q_func,
		std::vector<double>& stiff, std::vector<double>& vec) const{
	_THROW_NOT_IMP_;
}

void ASpatialApproximator::apply_bc_robin_to_stiff(
		int ibnd,
		std::function<double(Point)> alpha_func,
		std::function<double(Point)> beta_func,
		std::vector<double>& stiff, std::vector<double>& vec) const{
	_THROW_NOT_IMP_;
}

std::string ASpatialApproximator::_vtk_outgrid() const{
	return _grid->vtk_outgrid();
}

std::vector<double> ASpatialApproximator::_prepare_vtk_scalar(const double* x) const{
	int n = _vtk_use_cell_data() ? _grid->n_cells() : _grid->n_points();
	return std::vector<double>(x, x+n);
}

bool ASpatialApproximator::_vtk_use_cell_data() const{
	_THROW_NOT_IMP_;
}

void ASpatialApproximator::vtk_save_scalar(std::string fname, const std::vector<double>& x, std::string dataname) const{
	vtk_save_scalar(fname, {x.data()}, {dataname});
}

void ASpatialApproximator::vtk_save_scalar(std::string fname, std::vector<const double*> x, std::vector<std::string> datanames) const{
	std::ofstream ofs(fname);
	if (!ofs){
		throw std::runtime_error("Failed to open file " + fname + " for writing");
	}
	// write grid
	if (_cache.vtk_outgrid.empty()){
		_cache.vtk_outgrid = _vtk_outgrid();
	}
	ofs << _cache.vtk_outgrid;

	if (x.size() == 0) return;

	// write data
	for (size_t idata=0; idata<x.size(); ++idata){
		std::vector<double> outvec = _prepare_vtk_scalar(x[idata]);
		if (idata == 0){
			if (_vtk_use_cell_data()){
				ofs << "CELL_DATA " << outvec.size() << std::endl;
			} else {
				ofs << "POINT_DATA " << outvec.size() << std::endl;
			}
		}
		ofs << "SCALARS " << datanames[idata] << " double 1" << std::endl;
		ofs << "LOOKUP_TABLE default" << std::endl;
		for (auto v: outvec) ofs << v << std::endl;
	}
}
