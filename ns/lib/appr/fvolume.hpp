#ifndef FVOLUME_HPP
#define FVOLUME_HPP

#include "common.hpp"
#include "grid/grid.hpp"

struct FFace{
	const int grid_index;
	const int negative_cell;
	const int positive_cell;
	const double area;
	const Point center;
	const Point normal;

	static FFace build(int iface, const std::vector<Point>& coo, int neg_cell, int pos_cell);
private:
	FFace(int iface, int neg_cell, int pos_cell, double area, Point center, Point normal);
};

struct FVolume{
	const int grid_index;
	const double volume;
	const Point center;

	static FVolume build(int icell, const std::vector<Point>& coo, CellCode code);
	static FVolume build_boundary(int icell, const std::vector<Point>& coo);
private:
	FVolume(int icell, double volume, Point center);
};

#endif
