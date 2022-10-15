#include "appr/fvm_approximator.hpp"

std::shared_ptr<FvmApproximator> FvmApproximator::build(std::shared_ptr<Grid> grid){
	return std::shared_ptr<FvmApproximator>(new FvmApproximator(grid));
}

FvmApproximator::FvmApproximator(std::shared_ptr<Grid> grid): ASpatialApproximator(grid){
	// internal faces
	for (int iface=0; iface<grid->n_faces(); ++iface){
		if (grid->is_boundary_face(iface)){
			continue;
		}
		Grid::FaceCellEntry fc = grid->tab_face_cell(iface);
		std::vector<Point> points;
		for (int ipoint: grid->tab_face_point(iface)){
			points.push_back(grid->point(ipoint));
		}
		_internal_faces.push_back(FFace::build(iface, points, fc.left_cell, fc.right_cell));
	}

	// boundary faces
	for (int btype: grid->btypes()){
		const GridBoundary& bnd = grid->boundary(btype);
		std::vector<FFace>& bfaces = _boundary_faces[btype];

		const std::vector<int>& ifaces = bnd.face_indices();
		const std::vector<int>& icells = bnd.cell_indices();
		for (int iface_loc=0; iface_loc<(int)ifaces.size(); ++iface_loc){
			std::vector<int> ipoints = bnd.tab_face_point_positive(iface_loc);
			std::vector<Point> points;
			for (int ipoint: grid->tab_face_point(ifaces[iface_loc])){
				points.push_back(grid->point(ipoint));
			}
			bfaces.push_back(FFace::build(ifaces[iface_loc], points, icells[iface_loc], -1));
		}
	}
	
	// volumes
	for (int icell=0; icell<grid->n_cells(); ++icell){
		std::vector<Point> points;
		for (int ipoint: grid->tab_cell_point(icell)){
			points.push_back(grid->point(ipoint));
		}
		_volumes.push_back(FVolume::build(icell, points, _grid->cell_code(icell)));
	}
}

std::vector<double> FvmApproximator::_build_stiff() const{
	const CsrStencil& sten = stencil();
	std::vector<double> ret(sten.n_nonzero(), 0);

	for (auto& f: _internal_faces){
		int i = f.negative_cell;
		int j = f.positive_cell;

		Vector vec_ij = _volumes[i].center - _volumes[j].center;
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
	for (int i=0; i<grid().n_cells(); ++i){
		stencil_set[i].insert(i);
	}

	// non-diagonals
	for (int iface=0; iface<grid().n_faces(); ++iface){
		const Grid::FaceCellEntry& fc = grid().tab_face_cell(iface);
		if (fc.left_cell >=0 && fc.right_cell >=0){
			stencil_set[fc.left_cell].insert(fc.right_cell);
			stencil_set[fc.right_cell].insert(fc.left_cell);
		}
	}

	return CsrStencil::build(stencil_set);
}

std::vector<double> FvmApproximator::approximate(std::function<double(Point)> func) const{
	std::vector<double> ret;
	
	for (auto& vol: _volumes){
		ret.push_back(func(vol.center));
	}

	return ret;
}

void FvmApproximator::apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff) const{
	for (const FFace& face: _boundary_faces.find(bnd)->second){
		int ivol = face.negative_cell;
		Vector v = face.center - _volumes[ivol].center;
		double val = face.area/vector_len(v)*vector_cos(v, face.normal);
		stiff[stencil().addr()[ivol]] += val;
	}
}

void FvmApproximator::apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> func, std::vector<double>& vec) const{
	for (const FFace& face: _boundary_faces.find(bnd)->second){
		int ivol = face.negative_cell;
		Vector v = face.center - _volumes[ivol].center;
		double val = face.area/vector_len(v)*vector_cos(v, face.normal);
		vec[ivol] += val * func(face.center);
	}
}

bool FvmApproximator::_vtk_use_cell_data() const{
	return true;
}
