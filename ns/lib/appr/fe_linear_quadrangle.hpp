#ifndef FE_LINEAR_QUADRANGLE_HPP
#define FE_LINEAR_QUADRANGLE_HPP

#include "common.hpp"
#include "appr/fe_element.hpp"

class FeLinearQuadrangle: public ALagrangeInternalElement{
public:
	FeLinearQuadrangle(const std::vector<int>& indices, const std::vector<Point>& points);
	int n_bases() const override { return 4; };
};


#endif
