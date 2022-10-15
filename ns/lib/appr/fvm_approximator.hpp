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

	void apply_bc_dirichlet_to_stiff_mat(int bnd, std::vector<double>& stiff) const override;
	void apply_bc_dirichlet_to_stiff_vec(int bnd, std::function<double(Point)> func, std::vector<double>& vec) const override;

	static std::shared_ptr<FvmApproximator> build(std::shared_ptr<Grid> grid);
protected:
	std::vector<double> _build_stiff() const override;
	std::vector<double> _build_mass() const override;
	CsrStencil _build_stencil() const override;
private:
	FvmApproximator(std::shared_ptr<Grid> grid);
	std::vector<FVolume> _volumes;
	std::vector<FFace> _internal_faces;
	std::map<int, std::vector<FFace>> _boundary_faces;

	struct Cache{
	};
	mutable Cache _cache;

	bool _vtk_use_cell_data() const override;
};

#endif
