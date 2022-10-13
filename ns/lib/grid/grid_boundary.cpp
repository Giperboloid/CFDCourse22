#include "grid/grid_boundary.hpp"
#include "grid/grid.hpp"

int GridBoundary::n_faces() const{
	return (int)_faces.size();
}

std::vector<int> GridBoundary::tab_faces() const{
	return _faces;
}

std::vector<int> GridBoundary::tab_face_point_positive(int iface) const{
	int gface = _faces[iface];
	std::vector<int> ret = _grid->tab_face_point(gface);
	if (_grid->tab_face_cell(gface).left_cell == -1){
		std::reverse(ret.begin(), ret.end());
	}
	return ret;
}

std::vector<int> GridBoundary::tab_points() const{
	if (_cache.tab_points.size() == 0){
		std::set<int> pp;
		for (int iface: _faces){
			auto v = _grid->tab_face_point(iface);
			pp.insert(v[0]);
			pp.insert(v[1]);
		}
		_cache.tab_points.resize(pp.size());
		std::copy(pp.begin(), pp.end(), _cache.tab_points.begin());
	}
	return _cache.tab_points;
}
