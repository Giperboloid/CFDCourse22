#include <set>
#include "appr/fem_approximator.hpp"

void AFemApproximator::_initialize(){
	// internal cells
	for (int icell=0; icell<_grid->n_cells(); ++icell){
		std::vector<int> ipoints = _grid->tab_cell_point(icell);
		std::vector<Point> points;
		for (int ip: ipoints){
			points.push_back(_grid->point(ip));
		}
		std::shared_ptr<AInternalElement> el = _build_internal_element(ipoints, points, _grid->cell_code(icell));
		_elements.push_back(el);
	}

	// boundary cells
	for (int b: _grid->btypes()){
		const GridBoundary& bnd = _grid->boundary(b);
		for (int iface=0; iface<bnd.n_faces(); ++iface){
			std::vector<int> ipoints = bnd.tab_face_point_positive(iface);
			std::vector<Point> points;
			for (int ip: ipoints){
				points.push_back(_grid->point(ip));
			}
			std::shared_ptr<AFaceElement> el = _build_boundary_element(ipoints, points);
			_boundary_elements[b].push_back(el);
		}
	}

	// initialize internal elements
	const CsrStencil& sten = stencil();
	for (auto& elem: _elements) elem->_initialize_stencil_addresses(sten);

	// initialize boundary elements
	for (auto& bit: _boundary_elements)
	for (auto& elem: bit.second){
		elem->_initialize_stencil_addresses(stencil());
	}

	std::cout << "FEM initialized with " << n_bases() << " basis functions" << std::endl;
}

CsrStencil AFemApproximator::_build_stencil() const{
	// stencil
	std::vector<std::set<int>> data(n_bases());
	for (auto elem: _elements){
		for (int irow: elem->global_connect())
		for (int icol: elem->global_connect()){
			if (irow != icol){
				data[irow].insert(icol);
			}
		}
	}
	std::vector<int> addr {0};
	std::vector<int> cols;

	for (int irow = 0; irow < (int)data.size(); ++irow){
		cols.push_back(irow);
		for (int icol: data[irow]){
			cols.push_back(icol);
		}
		addr.push_back(addr.back() + data[irow].size() + 1);
	}

	CsrStencil ret;
	ret.set_data(std::move(addr), std::move(cols));

	return ret;
}

std::vector<double> AFemApproximator::_build_stiff() const{
	std::vector<double> ret(stencil().n_nonzero(), 0);
	for (auto el: _elements){
		el->add_to_global_mat(el->stiff(), ret);
	}
	return ret;
}

std::vector<double> AFemApproximator::_build_mass() const{
	std::vector<double> ret(stencil().n_nonzero(), 0);
	for (auto el: _elements){
		el->add_to_global_mat(el->mass(), ret);
	}
	return ret;
}
std::vector<double> AFemApproximator::_build_lumped_mass() const{
	_THROW_NOT_IMP_;
}

std::vector<Point> AFemApproximator::element_centered_gradient(const std::vector<double>& fun) const{
	std::vector<Point> ret;
	for (auto el: _elements){
		ret.push_back(el->element_centered_gradient(fun));
	}
	return ret;
}

std::vector<double> AFemApproximator::transport_v_elem(
		const std::vector<double>& vx,
		const std::vector<double>& vy,
		const std::vector<double>& vz) const {
	std::vector<double> ret(stencil().n_nonzero(), 0);
	bool need_x = !vx.empty();
	bool need_y = !vy.empty();
	bool need_z = !vz.empty();

	std::vector<double> grad_x, grad_y, grad_z;

	for (int ielem=0; ielem<n_elements(); ++ielem){
		const AInternalElement& elem = *_elements[ielem];
		if (need_x) grad_x = elem.grad_x();
		if (need_y) grad_y = elem.grad_y();
		if (need_z) grad_z = elem.grad_z();

		int k = 0;
		for (int irow: elem.global_connect())
		for (int icol: elem.global_connect()){
			int addr = stencil().addr_index(irow, icol);
			if (need_x) ret[addr] -= vx[ielem] * grad_x[k];
			if (need_y) ret[addr] -= vy[ielem] * grad_y[k];
			if (need_z) ret[addr] -= vz[ielem] * grad_z[k];
			++k;
		}
	}
	return ret;
}

const std::vector<int>& AFemApproximator::boundary_bases(int btype) const{
	if (_cache.boundary_bases.empty()){
		for (const auto& it: _boundary_elements){
			std::set<int> bv;
			for (const auto& belem: it.second){
				for (int iv: belem->global_connect()){
					bv.insert(iv);
				}
			}
			_cache.boundary_bases[it.first] = std::vector<int>(bv.begin(), bv.end());
		}
	}
	return _cache.boundary_bases[btype];
}
