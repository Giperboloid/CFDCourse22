#ifndef FE_LINEAR_SEGMENT_HPP
#define FE_LINEAR_SEGMENT_HPP

#include "common.hpp"
#include "appr/fe_element.hpp"

class FeLinearSegment: public ALagrangeFaceElement{
public:
	FeLinearSegment(const std::vector<int>& indices, const std::vector<Point>& points);

	int n_bases() const override { return 2; };
};



#endif
