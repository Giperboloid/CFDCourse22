#include <iostream>
#include <map>
#include <vector>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>

// ====================== CSRMatrix
struct CSRMatrix{
	std::vector<size_t> ptr;
	std::vector<size_t> col;
	std::vector<double> val;
};

// ======================= SparseMatrix
class SparseMatrix{
public:
	typedef std::map<size_t, double> row_t;

	// Constructs empty NxN matrix
	SparseMatrix(size_t N);

	// Returns number of rows
	size_t dim() const;
	// Returns number of non-zero entries
	size_t n_non_zero() const;

	// Returns value at (irw, icol) position
	double value(size_t irow, size_t icol) const;
	// Returns read-only irow-th row as {column: value} map.
	const row_t& row(size_t irow) const;

	// Clears irow-th row and sets unity to its diagonal position
	void unit_row(size_t irow);

	// M[irow, icol] += val
	void add(size_t irow, size_t icol, double val);

	// M[irow, icol] = 0
	void rem(size_t irow, size_t icol);

	// solves M*x = rhs using amgcl
	// returns <number of iterations, final tolerance> tuple
	std::tuple<size_t, double> solve(const std::vector<double>& rhs, std::vector<double>& x) const;

	// Performs matrix-vector multiplication in form:
	//     ret += coef * M*x
	void mult_vec(const std::vector<double>& x, double coef, std::vector<double>& ret) const;

	// Converts matrix to csr format
	CSRMatrix as_csr() const;
	// Converts matrix to dense format
	std::vector<double> as_dense() const;

	// Sets maximum allowed iterations and minimum allowed tolerance for the amgcl solver
	void set_solver_params(size_t maxit, double tolerance);
private:
	const size_t N;
	std::vector<row_t> data;
	size_t solver_maximum_iterations = 1000;
	double solver_tolerance = 1e-12;
};

// ======================== SparseMatrix implementation
std::ostream& operator<<(std::ostream& os, const SparseMatrix& M){
	if (M.dim() <= 12){
		// use dense output for small matrices
		std::vector<double> dt = M.as_dense();
		double* it = dt.data();
		for (size_t irow=0; irow < M.dim(); ++irow){
			for (size_t icol=0; icol < M.dim(); ++icol){
				os << std::setw(10) << *it++ << " ";
			}
			os << std::endl;
		}
	} else {
		// prints only non-zero entries for large matrices
		for (size_t irow=0; irow<M.dim(); ++irow){
			os << "--- ROW " << irow << std::endl;
			os << "/";
			for (auto& it: M.row(irow)){
				os << " " << it.first << ": " << it.second << " /";
			}
			os << std::endl;
		}
	}
	return os;
}

SparseMatrix::SparseMatrix(size_t N): N(N), data(N){ }

size_t SparseMatrix::dim() const{
	return N;
}
const SparseMatrix::row_t& SparseMatrix::row(size_t irow) const{
	return data[irow];
}

double SparseMatrix::value(size_t irow, size_t icol) const{
	auto fnd = data[irow].find(icol);
	if (fnd == data[irow].end()){
		return 0;
	} else {
		return fnd->second;
	}
}

size_t SparseMatrix::n_non_zero() const{
	size_t ret = 0;
	for (auto& row: data) ret += row.size();
	return ret;
}
void SparseMatrix::unit_row(size_t irow){
	data[irow].clear();
	data[irow].emplace(irow, 1);
}

void SparseMatrix::add(size_t irow, size_t icol, double val){
	row_t& row = data[irow];
	auto r = row.emplace(icol, val);
	if (!r.second) r.first->second += val;
}

void SparseMatrix::rem(size_t irow, size_t icol){
	row_t& row = data[irow];
	auto r = row.find(icol);
	if (r != row.end()){
		row.erase(r);
	}
}

void SparseMatrix::mult_vec(const std::vector<double>& x, double coef, std::vector<double>& ret) const{
	ret.resize(dim(), 0);
	for (size_t irow=0; irow<dim(); ++irow){
		for(auto& it: data[irow]){
			ret[irow] += coef * it.second * x[it.first];
		}
	}
}

void SparseMatrix::set_solver_params(size_t maxit, double tolerance){
	solver_maximum_iterations = maxit;
	solver_tolerance = tolerance;
}

CSRMatrix SparseMatrix::as_csr() const{
	CSRMatrix ret;
	size_t nn = n_non_zero();

	ret.ptr.resize(dim() + 1, 0);
	ret.col.resize(nn);
	ret.val.resize(nn);

	size_t k = 0;
	for (size_t irow=0; irow<dim(); ++irow){
		for (auto& it: data[irow]){
			ret.col[k] = it.first;
			ret.val[k] = it.second;
			++k;
		}
		ret.ptr[irow+1] = k;
	}

	return ret;
}

std::vector<double> SparseMatrix::as_dense() const{
	std::vector<double> ret(N*N, 0);
	for (size_t irow = 0; irow < N; ++irow){
		const row_t& row = data[irow];
		for (auto& it: row){
			size_t gind = N*irow + it.first;
			ret[gind] = it.second;
		}
	}
	return ret;
}
std::tuple<size_t, double> SparseMatrix::solve(const std::vector<double>& rhs, std::vector<double>& x) const{
	CSRMatrix csr = as_csr();

	amgcl::backend::crs<double, size_t> amgcl_matrix;
	amgcl_matrix.own_data = false;
	amgcl_matrix.nrows = amgcl_matrix.ncols = dim();
	amgcl_matrix.nnz = n_non_zero();
	amgcl_matrix.ptr = csr.ptr.data();
	amgcl_matrix.col = csr.col.data();
	amgcl_matrix.val = csr.val.data();

	boost::property_tree::ptree prm;
	prm.put("solver.type", "fgmres");
	prm.put("solver.tol", solver_tolerance);
	prm.put("solver.maxiter", solver_maximum_iterations);
	prm.put("precond.coarsening.type", "smoothed_aggregation");
	prm.put("precond.relax.type", "spai0");

	using backend_t = amgcl::backend::builtin<double>;
	using solver_t = amgcl::make_solver<
		amgcl::amg<
			backend_t,
			amgcl::runtime::coarsening::wrapper,
			amgcl::runtime::relaxation::wrapper
		>,
		amgcl::runtime::solver::wrapper<backend_t>
	>;
	solver_t slv(amgcl_matrix, prm);

	x.resize(dim());
	auto ret = slv(rhs, x);

	if (std::get<0>(ret) >= solver_maximum_iterations){
		std::ostringstream os;
		os << "WARNING: Sparse matrix solution failed to converge with tolerance: "
			<< std::scientific << std::setprecision(2) << std::get<1>(ret)
			<< std::endl;
		std::cout << os.str();
	}

	return ret;
}

int main(){
	size_t n = 10;

	// Fill M and rhs
	SparseMatrix M(n);
	std::vector<double> rhs(n, 0);

	for (size_t i=1; i<n-1; ++i){
		M.add(i, i, 2.0);
		M.add(i, i-1, -1.0);
		M.add(i, i+1, -1.0);
	}

	M.unit_row(0);
	rhs[0] = 1;
	M.unit_row(n-1);
	rhs[n-1] = 2;

	// Print M and rhs
	std::cout << "========== Matrix" << std::endl;
	std::cout << M << std::endl;
	std::cout << "========== Right hand side" << std::endl;
	for (size_t i=0; i<n; ++i){
		std::cout << rhs[i] << " ";
	}
	std::cout << std::endl;

	// Solve SLE
	std::vector<double> x(n);
	std::tuple<size_t, double> r = M.solve(rhs, x);

	// Print solution
	std::cout << "========== Solution" << std::endl;
	for (size_t i=0; i<n; ++i){
		std::cout << x[i] << " ";
	}
	std::cout << std::endl;

	// Print additional solution info
	std::cout << "========== Solution properties" << std::endl;
	std::cout << "iterations: " << std::get<0>(r) << std::endl;
	std::cout << "tolerance : " << std::get<1>(r) << std::endl;
}
