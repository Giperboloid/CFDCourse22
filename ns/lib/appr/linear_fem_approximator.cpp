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

bool LinearFemApproximator::_vtk_use_cell_data() const{
	return false;
}

std::map<int, std::vector<std::pair<int, Point>>> LinearFemApproximator::_build_boundary_bases() const{
	std::map<int, std::vector<std::pair<int, Point>>> ret;

	for (int ibnd: _grid->btypes()){
		auto fnd = _boundary_elements.find(ibnd);
		if (fnd == _boundary_elements.end()){
			throw std::runtime_error("Failed to find boundary " + std::to_string(ibnd));
		}
		std::vector<std::pair<int, Point>> p;

		for (std::shared_ptr<const AFaceElement> el: fnd->second){
			const ALagrangeFaceElement* face = static_cast<const ALagrangeFaceElement*>(el.get());
			const std::vector<int>& ind = face->global_connect();
			const std::vector<Point>& pts = face->coo();

			for (size_t i=0; i<ind.size(); ++i){
				p.push_back({ind[i], pts[i]});
			}
		}

		ret[ibnd] = p;
	}

	return ret;

	return ret;
}
