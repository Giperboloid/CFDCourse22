#ifndef FE_ELEMENT_HPP
#define FE_ELEMENT_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "geom.hpp"

// abstract finite element
class AElement{
	friend class AFemApproximator;
public:
	AElement(const std::vector<int>& connect): _connect(connect){}
	virtual ~AElement() = default;

	// number of local bases in the element
	virtual int n_bases() const = 0;

	// == assembly procedures
	// global basis indices of the element local basis functions
	const std::vector<int>& global_connect() const { return _connect; }
	// add dense local matrix to the global sparse fem matrix
	void add_to_global_mat(const std::vector<double>& lmat,
	                       std::vector<double>& mat_values);
	// add local vector to the global fem vector
	void add_to_global_vec(const std::vector<double>& lvec,
	                       std::vector<double>& vec);
	void clear_cache() const noexcept;
protected:
	std::vector<int> _stencil_addr;
	std::vector<int> _connect;
	mutable struct Cache {
		std::vector<double> grad_x;
		std::vector<double> grad_y;
		std::vector<double> grad_z;

		bool is_grad_x{false};
		bool is_grad_y{false};
		bool is_grad_z{false};
	} _cache;
private:
	void _initialize_stencil_addresses(const CsrStencil& M);
};

// Internal finite element
class AInternalElement: public AElement{
public:
	AInternalElement(const std::vector<int>& basis_indices): AElement(basis_indices){};

	// === local matrices
	// [i, j] = grad(p_i) @ grad(p_j)
	virtual std::vector<double> stiff() const;
	// [i, j] = p_i*p_j
	virtual std::vector<double> mass() const;
	// [i, j] = dp_j/dx * p_i
	virtual const std::vector<double>& grad_x() const;
	// [i, j] = dp_j/dy * p_i
	virtual const std::vector<double>& grad_y() const;
	// [i, j] = dp_j/dz * p_i
	virtual const std::vector<double>& grad_z() const;

	// clear cached data

	// [df/dx, df/dy, df/dz] at the element center
	virtual Point element_centered_gradient(const std::vector<double>& fun) const;
};

// Face finite element
class AFaceElement: public AElement{
public:
	AFaceElement(const std::vector<int>& basis_indices): AElement(basis_indices){};
};

// Lagrangian internal finite element
class ALagrangeInternalElement: public AInternalElement{
public:
	ALagrangeInternalElement(const std::vector<int>& vert_indices, const std::vector<Point>& vert_coo)
		: AInternalElement(vert_indices), _coo(vert_coo){};
protected:
	std::vector<Point> _coo;
};


// Lagrangian face finite element
class ALagrangeFaceElement: public AFaceElement{
public:
	ALagrangeFaceElement(const std::vector<int>& vert_indices, const std::vector<Point>& vert_coo)
		: AFaceElement(vert_indices), _coo(vert_coo){};
protected:
	std::vector<Point> _coo;
};



#endif
