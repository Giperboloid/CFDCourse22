#ifndef GRID2_BOUNDARY_HPP
#define GRID2_BOUNDARY_HPP

#include "common.hpp"

class Grid;

class GridBoundary{
	friend class GridBuilder;
public:
	// number of faces in the boundary
	int n_faces() const;
	// list of faces of the boundary
	std::vector<int> tab_faces() const;
	// boundary face points from the boundary face index.
	// Returned face is directed along the boundary
	std::vector<int> tab_face_point_positive(int ilocal_face) const;
	// list of boundary points
	std::vector<int> tab_points() const;
private:
	struct Cache{
		std::vector<int> tab_points;
	};
	mutable Cache _cache;

	GridBoundary(const Grid* grid, std::vector<int>&& faces): _grid(grid), _faces(std::move(faces)) {}

	const Grid* _grid;
	std::vector<int> _faces;
};


#endif
