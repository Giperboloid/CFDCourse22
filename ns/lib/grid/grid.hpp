#ifndef GRID_HPP
#define GRID_HPP

#include "common.hpp"
#include "geom.hpp"
#include "grid/grid_boundary.hpp"

enum struct CellCode{
	SEGMENT = 3,
	POLYGON = 7,
	TETRAHEDRON = 10,
	HEXAHEDRON = 12,
	WEDGE = 13,
	PYRAMID = 14,
	PENTAPRISM = 15,
	HEXAPRISM = 16
};

class Grid{
	friend class GridBuilder;
public:
	// === geometrical dimension
	const int dim;

	// === sizes
	int n_points() const;
	int n_faces() const;
	int n_cells() const;
	
	// === geometry
	Point point(int point_index) const;
	std::vector<int> btypes() const;
	const GridBoundary& boundary(int ibnd) const;
	bool is_boundary_face(int iface) const;
	CellCode cell_code(int icell) const;
	
	// === connectivities
	std::vector<int> tab_cell_point(int icell) const;
	std::vector<int> tab_point_cell(int icell) const;
	std::vector<int> tab_face_point(int iface) const;
	std::vector<int> tab_point_face(int iface) const;

	struct CellFaceEntry{
		int face_index;
		int normal_direction;
	};
	std::vector<CellFaceEntry> tab_cell_face(int icell) const;

	struct FaceCellEntry{
		int right_cell = -1;
		int left_cell = -1;
	};
	FaceCellEntry tab_face_cell(int iface) const;

	// === savers
	std::string vtk_outgrid() const;
	void save_vtk(std::string outfile) const;
	std::string vtk_outgrid_faces() const;
	void save_vtk_faces(std::string outfile) const;
private:
	struct Cache{
		std::vector<std::vector<CellFaceEntry>> tab_cell_face;
		std::vector<std::vector<int>> tab_point_cell;
		std::vector<std::vector<int>> tab_point_face;
	};
	mutable Cache _cache;

	Grid(int dim): dim(dim){};

	// data
	std::vector<Point> _points;
	std::vector<std::vector<int>> _cells;
	std::vector<std::vector<int>> _faces;
	std::vector<FaceCellEntry> _tab_face_cell;

	std::map<int, GridBoundary> _boundaries;

	std::vector<int> _vtk_cell_codes;
};


#endif
