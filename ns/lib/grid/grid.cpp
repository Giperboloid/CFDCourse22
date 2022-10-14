#include <list>
#include <fstream>
#include <sstream>
#include <numeric>
#include "grid/grid.hpp"

namespace{
std::vector<std::vector<int>> invert_tab(const std::vector<std::vector<int>>& tab){
	// find maximum
	int mx = -1;
	for (const auto& v: tab){
		mx = std::max(mx, *std::max_element(v.begin(), v.end()));
	}

	// allocate return
	std::vector<std::vector<int>> ret(mx+1);

	// invert
	for (int irow=0; irow<(int)tab.size(); ++irow){
		for (int icol: tab[irow]){
			ret[icol].push_back(irow);
		}
	}

	// remove duplicates
	for (auto& v: ret){
		std::sort(v.begin(), v.end());
		v.resize(std::unique(v.begin(), v.end()) - v.begin());
	}

	return ret;
}
};

int Grid::n_points() const{
	return _points.size();
}

int Grid::n_cells() const{
	return _cells.size();
}

int Grid::n_faces() const{
	return _faces.size();
}

std::string Grid::vtk_outgrid() const{
	std::ostringstream ofs;

	// header
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "cfdlib" << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;

	// points
	ofs << "POINTS " << _points.size() << " double" << std::endl;
	for (size_t i=0; i < _points.size(); ++i){
		ofs << _points[i].x << " " << _points[i].y << " " << _points[i].z << std::endl;
	}

	// cells
	int ntot = std::accumulate(_cells.begin(), _cells.end(), 0,
	                           [](int m, const std::vector<int>& v){return m + v.size(); })
	           + _cells.size();
	ofs << "CELLS " << n_cells() << " " << ntot << std::endl;
	for (int icell=0; icell<n_cells(); ++icell){
		ofs << _cells[icell].size();
		for (size_t ipoint=0; ipoint<_cells[icell].size(); ++ipoint){
			ofs << " " << _cells[icell][ipoint];
		}
		ofs << std::endl;
	}

	// celltypes
	ofs << "CELL_TYPES " << n_cells() << std::endl;
	for (int icell=0; icell<n_cells(); ++icell){
		ofs << _vtk_cell_codes[icell] << std::endl;
	}

	return ofs.str();
}

std::string Grid::vtk_outgrid_faces() const{
	std::ostringstream ofs;

	// header
	ofs << "# vtk DataFile Version 2.0" << std::endl;
	ofs << "cfdlib" << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;

	// points
	ofs << "POINTS " << _points.size() << " double" << std::endl;
	for (size_t i=0; i < _points.size(); ++i){
		ofs << _points[i].x << " " << _points[i].y << " " << _points[i].z << std::endl;
	}

	// cells
	int ntot = std::accumulate(_faces.begin(), _faces.end(), 0,
	                           [](int m, const std::vector<int>& v){return m + v.size(); })
	           + _faces.size();
	ofs << "CELLS " << n_faces() << " " << ntot << std::endl;
	for (int iface=0; iface < n_faces(); ++iface){
		ofs << _faces[iface].size();
		for (size_t ipoint=0; ipoint<_faces[iface].size(); ++ipoint){
			ofs << " " << _faces[iface][ipoint];
		}
		ofs << std::endl;
	}

	// celltypes
	ofs << "CELL_TYPES " << n_faces() << std::endl;
	int code = (dim == 2) ? 3 : 7;
	for (int iface=0; iface<n_faces(); ++iface){
		ofs << code << std::endl;
	}
	return ofs.str();
}

void Grid::save_vtk(std::string outfile) const{
	std::cout << "Writing vtk grid to " << outfile << std::endl;

	std::ofstream ofs(outfile);
	if (!ofs){
		throw std::runtime_error("failed to open " + outfile);
	}

	ofs << vtk_outgrid();
}

void Grid::save_vtk_faces(std::string outfile) const{
	std::cout << "Writing vtk grid faces to " << outfile << std::endl;

	std::ofstream ofs(outfile);
	if (!ofs){
		throw std::runtime_error("failed to open " + outfile);
	}

	ofs << vtk_outgrid_faces();

	// boundary type
	std::vector<int> bt(n_faces(), -1);
	for (int iface=0; iface<n_faces(); ++iface){
		if (is_boundary_face(iface)) bt[iface] = 0;
	}
	for (auto& bit: _boundaries){
		for (int iface: bit.second.tab_faces()){
			bt[iface] = bit.first;
		}
	}
	ofs << "CELL_DATA " << n_faces() << std::endl;
	ofs << "SCALARS btype int 1" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	for (int& v: bt){
		ofs << v << std::endl;
	}

}

Point Grid::point(int point_index) const{
	return _points[point_index];
}

std::vector<int> Grid::btypes() const{
	std::vector<int> ret;
	for (auto& it: _boundaries){ ret.push_back(it.first); }
	return ret;
}

const GridBoundary& Grid::boundary(int ibnd) const{
	auto fnd = _boundaries.find(ibnd);
	if (fnd == _boundaries.end()){
		throw std::runtime_error("Failed to find boundary with id = " + std::to_string(ibnd));
	}
	return fnd->second;
}

bool Grid::is_boundary_face(int iface) const{
	const auto& fc = _tab_face_cell[iface];
	return fc.left_cell < 0 || fc.right_cell < 0;
}

CellCode Grid::cell_code(int icell) const{
	return (CellCode)_vtk_cell_codes[icell];
}

std::vector<int> Grid::tab_cell_point(int icell) const{
	return _cells[icell];
}

std::vector<int> Grid::tab_face_point(int iface) const{
	return _faces[iface];
}

std::vector<int> Grid::tab_point_cell(int ipoint) const{
	if (_cache.tab_point_cell.empty()){
		_cache.tab_point_cell = invert_tab(_cells);
	}
	return _cache.tab_point_cell[ipoint];
}

std::vector<int> Grid::tab_point_face(int ipoint) const{
	if (_cache.tab_point_face.empty()){
		_cache.tab_point_face = invert_tab(_faces);
	}
	return _cache.tab_point_face[ipoint];
}

Grid::FaceCellEntry Grid::tab_face_cell(int iface) const{
	return _tab_face_cell[iface];
}

std::vector<Grid::CellFaceEntry> Grid::tab_cell_face(int icell) const{
	if (_cache.tab_cell_face.empty()){
		auto& tab = _cache.tab_cell_face;
		tab.resize(n_cells());

		for (int iface=0; iface<n_faces(); ++iface){
			const FaceCellEntry& fc = _tab_face_cell[iface];
			if (fc.left_cell >= 0){
				CellFaceEntry e;
				e.face_index = iface;
				e.normal_direction = 1;
				tab[fc.left_cell].push_back(e);
			}
			if (fc.right_cell >= 0){
				CellFaceEntry e;
				e.face_index = iface;
				e.normal_direction = -1;
				tab[fc.right_cell].push_back(e);
			}
		}
	}
	return _cache.tab_cell_face[icell];
}
