#ifndef GRID_HPP
#define GRID_HPP

#include "common.hpp"
#include "geom.hpp"
#include "grid/grid_boundary.hpp"

// valid cell types
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
	// point at index
	Point point(int point_index) const;
	// passed boundary types
	std::vector<int> btypes() const;
	// returns boundary for the specified type
	const GridBoundary& boundary(int ibnd) const;
	// if face with given index has only one adjacent cell
	bool is_boundary_face(int iface) const;
	// type of the cell
	CellCode cell_code(int icell) const;
	
	// === connectivities
	// indices of points adjacent of the given cell
	std::vector<int> tab_cell_point(int icell) const;
	// indices of cells adjacent to the given point
	std::vector<int> tab_point_cell(int ipoint) const;
	// indices of points adjacent to the given face
	std::vector<int> tab_face_point(int iface) const;
	// indices of faces adjacent to the given point
	std::vector<int> tab_point_face(int ipoint) const;

	// cell-face connection info
	struct CellFaceEntry{
		// index of the face
		int face_index;
		//  1 - if face normal is directed outside the cell
		// -1 - otherwise
		int normal_direction;
	};
	// faces adjacent to the given cell
	std::vector<CellFaceEntry> tab_cell_face(int icell) const;

	// face-cell connection info
	struct FaceCellEntry{
		// face normal direction:
		//   1d - along x-axis
		//   2d - to the right from the [point0, point1] vector
		//   3d - direction of the vector product [point0, point1]X[point1, point2] vectors
		//
		//   where point0, point1, ... - consecutive face points

		// index of the cell along the face normal direction
		// -1 if no cell (boundary face)
		int right_cell = -1;
		// index of the cell against the face normal direction
		// -1 if no cell (boundary face)
		int left_cell = -1;
	};
	// cells adjacent to the given face
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
	// [ cell->points ]
	std::vector<std::vector<int>> _cells;
	// [ face->points ]
	std::vector<std::vector<int>> _faces;
	// [ face->cells ]
	std::vector<FaceCellEntry> _tab_face_cell;

	// grid boundaries with the given btype
	std::map<int, GridBoundary> _boundaries;

	// cell types
	std::vector<int> _vtk_cell_codes;
};


#endif
