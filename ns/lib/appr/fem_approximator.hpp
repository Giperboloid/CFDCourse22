#ifndef FEM_APPROXIMATOR
#define FEM_APPROXIMATOR

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "appr/spatial_approximator.hpp"
#include "appr/fe_element.hpp"

class AFemApproximator: public ASpatialApproximator{
public:
	virtual ~AFemApproximator() = default;

	int n_elements() const { return (int)_elements.size(); }

	std::vector<Point> element_centered_gradient(const std::vector<double>& fun) const;
	const std::vector<int>& boundary_bases(int btype) const;

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


	std::vector<double> _build_stiff() const override;
	std::vector<double> _build_mass() const override;
	std::vector<double> _build_lumped_mass() const override;
	CsrStencil _build_stencil() const override;
	

	std::vector<std::shared_ptr<AInternalElement>> _elements;
	std::map<int, std::vector<std::shared_ptr<AFaceElement>>> _boundary_elements;
private:
	struct Cache{
		std::map<int, std::vector<int>> boundary_bases;
	};
	mutable Cache _cache;
};

#endif
