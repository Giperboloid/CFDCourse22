#ifndef FE_ELEMENT_HPP
#define FE_ELEMENT_HPP

#include "common.hpp"
#include "slae/csrmat.hpp"
#include "geom.hpp"

class AElement{
	friend class AFemApproximator;
public:
	AElement(const std::vector<int>& connect): _connect(connect){}
	virtual ~AElement() = default;
	virtual int n_bases() const = 0;

	// == assembly procedures
	const std::vector<int> global_connect() const { return _connect; }
	void add_to_global_mat(const std::vector<double>& lmat,
	                       std::vector<double>& mat_values);
	void add_to_global_vec(const std::vector<double>& lvec,
	                       std::vector<double>& vec);
protected:
	std::vector<int> _stencil_addr;
	std::vector<int> _connect;
private:
	void _initialize_stencil_addresses(const CsrStencil& M);
};

class AInternalElement: public AElement{
public:
	AInternalElement(const std::vector<int>& vert_indices): AElement(vert_indices){};

	virtual std::vector<double> stiff() const;
	virtual std::vector<double> mass() const;
	virtual std::vector<double> grad_x() const;
	virtual std::vector<double> grad_y() const;
	virtual std::vector<double> grad_z() const;

	virtual Point element_centered_gradient(const std::vector<double>& fun) const;
};

class AFaceElement: public AElement{
public:
	AFaceElement(const std::vector<int>& vert_indices): AElement(vert_indices){};
};

class ALagrangeInternalElement: public AInternalElement{
public:
	ALagrangeInternalElement(const std::vector<int>& vert_indices, const std::vector<Point>& vert_coo)
		: AInternalElement(vert_indices), _coo(vert_coo){};
protected:
	std::vector<Point> _coo;
};


class ALagrangeFaceElement: public AFaceElement{
public:
	ALagrangeFaceElement(const std::vector<int>& vert_indices, const std::vector<Point>& vert_coo)
		: AFaceElement(vert_indices), _coo(vert_coo){};
protected:
	std::vector<Point> _coo;
};



#endif
