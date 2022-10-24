#ifndef GRID2_BUILDER_HPP
#define GRID2_BUILDER_HPP

#include "common.hpp"
#include "grid/grid.hpp"

class GridBuilder{
public:
	// grid from the vtk grid file exported from the gmsh
	static std::shared_ptr<Grid> build_from_gmshvtk(std::string fn);
	static std::shared_ptr<Grid> build_regular1(double len, int n_points);
};

#endif
