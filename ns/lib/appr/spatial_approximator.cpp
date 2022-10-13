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
	if (_cache.stiff.size() == 0){
		_cache.stiff = _build_stiff();
	}
	return _cache.stiff;
}
const std::vector<double>& ASpatialApproximator::mass() const{
	if (_cache.mass.size() == 0){
		_cache.mass = _build_mass();
	}
	return _cache.mass;
}
const std::vector<double>& ASpatialApproximator::lumped_mass() const{
	if (_cache.lumped_mass.size() == 0){
		_cache.lumped_mass = _build_lumped_mass();
	}
	return _cache.lumped_mass;
}

std::string ASpatialApproximator::_vtk_outgrid() const{
	return _grid->vtk_outgrid();
}

std::vector<double> ASpatialApproximator::_prepare_vtk_scalar(const double* x) const{
	return std::vector<double>(x, x+n_bases());
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
