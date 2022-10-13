#include "appr/linear_fem_approximator.hpp"
#include "appr/fe_linear_triangle.hpp"
#include "appr/fe_linear_quadrangle.hpp"
#include "appr/fe_linear_segment.hpp"

std::shared_ptr<LinearFemApproximator> LinearFemApproximator::build(std::shared_ptr<Grid> grid){
	std::shared_ptr<LinearFemApproximator> ret(new LinearFemApproximator(grid));
	ret->_initialize();

	return ret;
}

LinearFemApproximator::LinearFemApproximator(std::shared_ptr<Grid> grid): AFemApproximator(grid){
}

int LinearFemApproximator::n_bases() const{
	return _grid->n_points();
};

std::vector<double> LinearFemApproximator::approximate(std::function<double(Point)> func) const{
	std::vector<double> ret(n_bases());

	for (int iv=0; iv<n_bases(); ++iv){
		ret[iv] = func(_grid->point(iv));
	}

	return ret;
}

std::shared_ptr<AInternalElement> LinearFemApproximator::_build_internal_element(
		const std::vector<int>& point_indices,
		const std::vector<Point>& points,
		CellCode code){

	if (code == CellCode::POLYGON){
		if (point_indices.size() == 3){
			return std::make_shared<FeLinearTriangle>(point_indices, points);
		} else if (point_indices.size() == 4){
			return std::make_shared<FeLinearQuadrangle>(point_indices, points);
		}
	}

	throw std::runtime_error("Linear 2d FEM grid should contain triangles and quadrangles only");
}

std::shared_ptr<AFaceElement> LinearFemApproximator::_build_boundary_element(
		const std::vector<int>& point_indices,
		const std::vector<Point>& points){

	if (point_indices.size() == 2){
		return std::make_shared<FeLinearSegment>(point_indices, points);
	} else {
		throw std::runtime_error("Linear 2d FEM grid should contain only linear boundary segments");
	}
}

void LinearFemApproximator::apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff) const{
	auto s = stencil();
	for (int ivert: _grid->boundary(bnd).tab_points()){
		s.set_unit_diagonal(ivert, stiff);
	}
}

void LinearFemApproximator::apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> func, std::vector<double>& vec) const{
	for (int ivert: _grid->boundary(bnd).tab_points()){
		vec[ivert] = func(_grid->point(ivert));
	}
}

bool LinearFemApproximator::_vtk_use_cell_data() const{
	return false;
}
