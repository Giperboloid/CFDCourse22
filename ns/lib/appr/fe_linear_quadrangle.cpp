#include "appr/fe_linear_quadrangle.hpp"

FeLinearQuadrangle::FeLinearQuadrangle(const std::vector<int>& indices, const std::vector<Point>& points)
		: ALagrangeInternalElement(indices, points){
}
