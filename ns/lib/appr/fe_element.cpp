#include "appr/fe_element.hpp"

void AElement::_initialize_stencil_addresses(const CsrStencil& M){
	_stencil_addr.resize(0);

	for (size_t irow: _connect)
	for (size_t icol: _connect){
		_stencil_addr.push_back(M.addr_index(irow, icol));
	}
}

void AElement::add_to_global_mat(const std::vector<double>& lmat,
                                 std::vector<double>& mat_values){
	for (size_t i=0; i<lmat.size(); ++i){
		mat_values[_stencil_addr[i]] += lmat[i];
	}
}

void AElement::add_to_global_vec(const std::vector<double>& lvec,
                                 std::vector<double>& vec){
	for (size_t i=0; i<_connect.size(); ++i){
		vec[_connect[i]] += lvec[i];
	}
}

void AElement::clear_cache() const noexcept {
	_cache.is_grad_x = false;
	_cache.is_grad_y = false;
	_cache.is_grad_z = false;

	_cache.grad_x.clear();
	_cache.grad_y.clear();
	_cache.grad_z.clear();
}

std::vector<double> AInternalElement::stiff() const{
	_THROW_NOT_IMP_;
}

std::vector<double> AInternalElement::mass() const{
	_THROW_NOT_IMP_;
}

const std::vector<double>& AInternalElement::grad_x() const{
	_THROW_NOT_IMP_;
}

const std::vector<double>& AInternalElement::grad_y() const{
	_THROW_NOT_IMP_;
}

const std::vector<double>& AInternalElement::grad_z() const{
	_THROW_NOT_IMP_;
}

Point AInternalElement::element_centered_gradient(const std::vector<double>& fun) const{
	_THROW_NOT_IMP_;
}
