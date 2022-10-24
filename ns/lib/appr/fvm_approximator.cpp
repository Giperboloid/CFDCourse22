#include "appr/fvm_approximator.hpp"

std::shared_ptr<FvmApproximator> FvmApproximator::build(std::shared_ptr<Grid> grid){
	return std::shared_ptr<FvmApproximator>(new FvmApproximator(grid));
}

FvmApproximator::FvmApproximator(std::shared_ptr<Grid> grid): ASpatialApproximator(grid){
	// build face->btype array
	std::vector<int> face_btype(grid->n_faces(), -1);
	for (int btype: grid->btypes()){
		for (int ifaces: grid->boundary(btype).face_indices()){
			face_btype[ifaces] = btype;
		}
	}

	// internal volumes
	for (int icell=0; icell<grid->n_cells(); ++icell){
		std::vector<Point> points;
		for (int ipoint: grid->tab_cell_point(icell)){
			points.push_back(grid->point(ipoint));
		}
		_volumes.push_back(FVolume::build(icell, points, _grid->cell_code(icell)));
	}

	// faces
	for (int iface=0; iface<grid->n_faces(); ++iface){
		Grid::FaceCellEntry fc = grid->tab_face_cell(iface);
		std::vector<Point> points;
		for (int ipoint: grid->tab_face_point(iface)){
			points.push_back(grid->point(ipoint));
		}

		// add boundary volumes
		if (fc.left_cell < 0 || fc.right_cell < 0){
			int btype = face_btype[iface];
			if (fc.left_cell < 0) {
				fc.left_cell = _volumes.size();
			} else {
				fc.right_cell = _volumes.size();
			}
			_boundary_volumes[btype].push_back(_volumes.size());
			_volumes.push_back(FVolume::build_boundary(_volumes.size(), points));
		}

		_faces.push_back(FFace::build(iface, points, fc.left_cell, fc.right_cell));
	}
}

std::vector<double> FvmApproximator::_build_stiff() const{
	const CsrStencil& sten = stencil();
	std::vector<double> ret(sten.n_nonzero(), 0);

	for (auto& f: _faces){
		int i = f.negative_cell;
		int j = f.positive_cell;

		Vector vec_ij = _volumes[j].center - _volumes[i].center;
		double cos_n_ij = vector_cos(f.normal, vec_ij);
		double h = vector_len(vec_ij);
		double value = f.area/h*cos_n_ij;

		// diagonals
		ret[sten.addr_index(i, i)] += value;
		ret[sten.addr_index(j, j)] += value;

		// non diagonals
		ret[sten.addr_index(i, j)] -= value;
		ret[sten.addr_index(j, i)] -= value;
	}

	return ret;
}

std::vector<double> FvmApproximator::_build_mass() const{
	const CsrStencil& sten = stencil();

	std::vector<double> ret(sten.n_nonzero(), 0);
	for (int ivol=0; ivol<n_bases(); ++ivol){
		int k = sten.addr()[ivol];
		ret[k] = _volumes[ivol].volume;
	}
	return ret;
}

CsrStencil FvmApproximator::_build_stencil() const{
	std::vector<std::set<int>> stencil_set(n_bases());

	// diagonals
	for (int i=0; i<n_bases(); ++i){
		stencil_set[i].insert(i);
	}

	// non-diagonals
	for (const FFace& face: _faces){
		stencil_set[face.positive_cell].insert(face.negative_cell);
		stencil_set[face.negative_cell].insert(face.positive_cell);
	}

	return CsrStencil::build(stencil_set);
}

std::map<int, std::vector<std::pair<int, Point>>> FvmApproximator::_build_boundary_bases() const{
	std::map<int, std::vector<std::pair<int, Point>>> ret;

	for (int ibnd: _grid->btypes()){
		auto fnd = _boundary_volumes.find(ibnd);
		if (fnd == _boundary_volumes.end()){
			throw std::runtime_error("Failed to find boundary " + std::to_string(ibnd));
		}
		std::vector<std::pair<int, Point>> p;

		for (int ivol: fnd->second){
			const FVolume& vol = _volumes[ivol];
			p.push_back({vol.grid_index, vol.center});
		}

		ret[ibnd] = p;
	}

	return ret;
}

std::vector<double> FvmApproximator::approximate(std::function<double(Point)> func) const{
	std::vector<double> ret;
	
	for (auto& vol: _volumes){
		ret.push_back(func(vol.center));
	}

	return ret;
}

bool FvmApproximator::_vtk_use_cell_data() const{
	return true;
}
