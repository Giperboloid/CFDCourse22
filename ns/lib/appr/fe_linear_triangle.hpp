#ifndef FE_LINEAR_TRIANGLE_HPP
#define FE_LINEAR_TRIANGLE_HPP

#include "common.hpp"
#include "appr/fe_element.hpp"


class FeLinearTriangle: public ALagrangeInternalElement{
public:
	FeLinearTriangle(const std::vector<int>& indices, const std::vector<Point>& points);

	int n_bases() const override { return 3; };

	std::vector<double> stiff() const override;
	std::vector<double> mass() const override;
	const std::vector<double>& grad_x() const override;
	const std::vector<double>& grad_y() const override;
	const std::vector<double>& grad_z() const override;

	Point element_centered_gradient(const std::vector<double>& fun) const override;
private:
	double _modj, _j21, _j22, _j11, _j12;
};

#endif
