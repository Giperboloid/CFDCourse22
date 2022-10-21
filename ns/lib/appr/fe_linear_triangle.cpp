#include "appr/fe_linear_triangle.hpp"

FeLinearTriangle::FeLinearTriangle(const std::vector<int>& indices, const std::vector<Point>& points)
		: ALagrangeInternalElement(indices, points){

	_j11 = points[1].x - points[0].x;
	_j12 = points[2].x - points[0].x;
	_j21 = points[1].y - points[0].y;
	_j22 = points[2].y - points[0].y;
	_modj = _j22*_j11-_j12*_j21;
}

std::vector<double> FeLinearTriangle::stiff() const{
	std::vector<double> ret(9);
	double j1222 = _j21 - _j22, j1121 = _j11 - _j12;

	ret[0] = ( j1222*j1222 + j1121*j1121)/2/_modj;
	ret[1] = ( _j22*j1222  + _j12*j1121 )/2/_modj;
	ret[2] = (-_j11*j1121  - _j21*j1222 )/2/_modj;
	ret[4] = ( _j22*_j22   + _j12*_j12  )/2/_modj;
	ret[5] = (-_j21*_j22   - _j12*_j11  )/2/_modj;
	ret[8] = ( _j11*_j11   + _j21*_j21  )/2/_modj;
	ret[3] = ret[1];
	ret[6] = ret[2];
	ret[7] = ret[5];

	return ret;
}

std::vector<double> FeLinearTriangle::mass() const{
	double v = _modj/24.0;
	return {2*v, v, v,
		v, 2*v, v,
		v, v, 2*v};
}

Point FeLinearTriangle::element_centered_gradient(const std::vector<double>& fun) const{
	double x = fun[_connect[0]]*(-_j22+_j21)/_modj
		+ fun[_connect[1]]*(_j22)/_modj
		+ fun[_connect[2]]*(-_j21)/_modj;

	double y = fun[_connect[0]]*(_j12-_j11)/_modj
		+ fun[_connect[1]]*(-_j12)/_modj
		+ fun[_connect[2]]*(_j11)/_modj;

	return {x, y, 0};
}

const std::vector<double>& FeLinearTriangle::grad_x() const{
	if (_cache.is_grad_x) {
		return _cache.grad_x;
	}
	_cache.grad_x.resize(9);
	double dx[3]={-_j22 + _j21, _j22, -_j21};

	int k=0;
	for (int i=0; i<3; ++i)
	for (int j=0; j<3; ++j) _cache.grad_x[k++] = dx[j] / 6.0;

	_cache.is_grad_x = true;

	return _cache.grad_x;
}

const std::vector<double>& FeLinearTriangle::grad_y() const{
	if (_cache.is_grad_y) {
		return _cache.grad_y;
	}
	_cache.grad_y.resize(9);
	double dx[3]={_j12 - _j11, -_j12, _j11};

	int k=0;
	for (int i=0; i<3; ++i)
	for (int j=0; j<3; ++j) _cache.grad_y[k++] = dx[j] / 6.0;

	_cache.is_grad_y = true;

	return _cache.grad_y;
}

const std::vector<double>& FeLinearTriangle::grad_z() const{
	_cache.grad_z.resize(9);
	return _cache.grad_z;
}
