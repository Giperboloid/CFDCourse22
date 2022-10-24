#ifndef FVM_APPROXIMATOR
#define FVM_APPROXIMATOR

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "appr/spatial_approximator.hpp"
#include "appr/fvolume.hpp"

// abstract finite element approximator
class FvmApproximator: public ASpatialApproximator{
public:
	virtual ~FvmApproximator() = default;

	int n_bases() const override { return (int)_volumes.size(); }
	
	std::vector<double> approximate(std::function<double(Point)> func) const override;

	static std::shared_ptr<FvmApproximator> build(std::shared_ptr<Grid> grid);
private:
	FvmApproximator(std::shared_ptr<Grid> grid);

	std::vector<FVolume> _volumes;
	std::map<int, std::vector<int>> _boundary_volumes;
	std::vector<FFace> _faces;

	std::vector<double> _build_stiff() const override;
	std::vector<double> _build_mass() const override;
	CsrStencil _build_stencil() const override;
	std::map<int, std::vector<std::pair<int, Point>>> _build_boundary_bases() const override;
	bool _vtk_use_cell_data() const override;
};

#endif
