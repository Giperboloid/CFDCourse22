#ifndef CSRMAT_HPP
#define CSRMAT_HPP

#include "common.hpp"

class CsrStencil{
public:
	void set_data(std::vector<int>&& addr, std::vector<int>&& cols);

	int n_nonzero() const;
	int n_rows() const;

	const std::vector<int>& addr() const { return _addr; }
	const std::vector<int>& cols() const { return _cols; }

	int addr_index(int irow, int jcol) const;
	void matvec(const std::vector<double>& mat, const std::vector<double>& vec, std::vector<double>& ret) const;
	void matvec_plus(double coeff, const std::vector<double>& mat, const std::vector<double>& vec, std::vector<double>& ret) const;
	double matvec_irow(int irow, const std::vector<double>& mat, const std::vector<double>& vec) const;
	void set_unit_diagonal(int irow, std::vector<double>& mat) const;
private:
	std::vector<int> _addr = {0};
	std::vector<int> _cols;
};


class CsrMatrix: public CsrStencil{
public:
	std::vector<double>& vals() { return _vals; }
	const std::vector<double>& vals() const { return _vals; }

	void matvec(const std::vector<double>& vec, std::vector<double>& ret) const;
	void matvec_plus(double coeff, const std::vector<double>& vec, std::vector<double>& ret) const;
	double matvec_irow(int irow, const std::vector<double>& vec) const;
	void set_unit_diagonal(int irow);
private:
	std::vector<double> _vals;
};


#endif
