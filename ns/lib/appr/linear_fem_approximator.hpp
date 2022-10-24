#ifndef LINEAR_FEM_APPROXIMATOR_HPP
#define LINEAR_FEM_APPROXIMATOR_HPP

#include "common.hpp"
#include "appr/fem_approximator.hpp"

class LinearFemApproximator: public AFemApproximator{
public:
	int n_bases() const override;

	std::vector<double> approximate(std::function<double(Point)> func) const override;

	static std::shared_ptr<LinearFemApproximator> build(std::shared_ptr<Grid> grid);
protected:
	LinearFemApproximator(std::shared_ptr<Grid> grid);

	std::shared_ptr<AInternalElement> _build_internal_element(
		const std::vector<int>& point_indices,
		const std::vector<Point>& points,
		CellCode code) override;

	std::shared_ptr<AFaceElement> _build_boundary_element(
		const std::vector<int>& point_indices,
		const std::vector<Point>& points) override;
private:
	bool _vtk_use_cell_data() const override;
	std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const override;
};



#endif
