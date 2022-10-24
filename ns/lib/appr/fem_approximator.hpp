#ifndef FEM_APPROXIMATOR
#define FEM_APPROXIMATOR

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "appr/spatial_approximator.hpp"
#include "appr/fe_element.hpp"

// abstract finite element approximator
class AFemApproximator: public ASpatialApproximator{
public:
	virtual ~AFemApproximator() = default;

	// number of finite elements
	int n_elements() const { return (int)_elements.size(); }

	// gradient of the fem approximated function in the element centers
	std::vector<Point> element_centered_gradient(const std::vector<double>& fun) const;

	// => matrix of the operator: vx*grad_x + vy*grad_y + vz*grad_z
	// where v is given as element centered data.
	//
	// If any of velocity component is zero, pass empty vector.
	// For example in 2d case: vz={}
	std::vector<double> transport_v_elem(
			const std::vector<double>& vx,
			const std::vector<double>& vy,
			const std::vector<double>& vz) const;
protected:
	AFemApproximator(std::shared_ptr<Grid> grid): ASpatialApproximator(grid){}
	void _initialize();

	virtual std::shared_ptr<AInternalElement> _build_internal_element(
		const std::vector<int>& point_indices,
		const std::vector<Point>& points,
		CellCode code) = 0;

	virtual std::shared_ptr<AFaceElement> _build_boundary_element(
		const std::vector<int>& point_indices,
		const std::vector<Point>& points) = 0;
	

	std::vector<std::shared_ptr<AInternalElement>> _elements;
	std::map<int, std::vector<std::shared_ptr<AFaceElement>>> _boundary_elements;
private:
	struct Cache{
	};
	mutable Cache _cache;

	std::vector<double> _build_stiff() const override;
	std::vector<double> _build_mass() const override;
	CsrStencil _build_stencil() const override;
};

#endif
