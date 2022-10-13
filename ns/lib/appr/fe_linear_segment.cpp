#include "appr/fe_linear_segment.hpp"

FeLinearSegment::FeLinearSegment(const std::vector<int>& indices, const std::vector<Point>& points)
		: ALagrangeFaceElement(indices, points){
}
