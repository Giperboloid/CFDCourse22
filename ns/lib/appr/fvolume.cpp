#include "appr/fvolume.hpp"

namespace{

std::pair<Point, double> face_normal_area(const std::vector<Point>& coo){
	if (coo.size() == 1){
		// point
		return { {1.0, 0.0, 0.0}, 1.0 };
	} else if (coo.size() == 2){
		// segment
		Vector v = coo[1] - coo[0];
		Vector ret = {v.y, -v.x, 0};
		double ret_len = vector_len(ret);
		return { ret/ret_len, ret_len };
	} else {
		// polygon
		Point ret;
		for (size_t i=1; i<coo.size()-1; ++i){
			Vector v1 = coo[i] - coo[0];
			Vector v2 = coo[i+1] - coo[0];
			ret += vector_cross(v1, v2);
		}
		double ret_len = vector_len(ret);
		return { ret/ret_len, ret_len/2 };
	}
}

Point aver_point(const std::vector<Point>& coo){
	Point ret;
	for (auto& p: coo) ret += p;
	ret /= coo.size();
	return ret;
}

double cell_volume(const std::vector<Point>& coo, CellCode code){
	switch (code){
		case CellCode::SEGMENT:
			return std::abs(coo[1].x - coo[0].x);
		case CellCode::POLYGON:
			return face_normal_area(coo).second;
		case CellCode::TETRAHEDRON:
			return vector_dot_triple(
				coo[1] - coo[0],
				coo[2] - coo[0],
				coo[3] - coo[0])/6;
		default:
			_THROW_NOT_IMP_;
	}
}

}

FFace::FFace(int iface, int neg_cell, int pos_cell, double area, Point center, Point normal):
	grid_index(iface),
	negative_cell(neg_cell),
	positive_cell(pos_cell),
	area(area),
	center(center),
	normal(normal) { }

FFace FFace::build(int iface, const std::vector<Point>& coo, int neg_cell, int pos_cell){
	Point normal;
	double area;
	std::tie(normal, area) = face_normal_area(coo);
	Point center = aver_point(coo);
	return FFace(iface, neg_cell, pos_cell, area, center, normal);
}

FVolume::FVolume(int icell, double volume, Point center):
	grid_index(icell),
	volume(volume),
	center(center){ }

FVolume FVolume::build(int icell, const std::vector<Point>& coo, CellCode code){
	Point center = aver_point(coo);
	double volume = cell_volume(coo, code);
	return FVolume(icell, volume, center);
}

FVolume FVolume::build_boundary(int icell, const std::vector<Point>& coo){
	return FVolume(icell, 0, aver_point(coo));
}
