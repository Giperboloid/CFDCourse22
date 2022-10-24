#include "dbg/printer.hpp"
#include <sstream>
#include <iomanip>

void dbg::print_matrow(int irow, const CsrStencil& stencil, const std::vector<double>& mat){
	std::cout << "==== matrix row = " << irow << " =========== " << std::endl;
	std::ostringstream ocol, oval;

	ocol << "cols: ";
	oval << "vals: ";

	for (int a=stencil.addr()[irow]; a < stencil.addr()[irow+1]; ++a){
		ocol << std::setw(14) << stencil.cols()[a] << " ";
		oval << std::setw(14) << std::setprecision(10) << mat[a] << " ";
	}

	std::cout << ocol.str() << std::endl;
	std::cout << oval.str() << std::endl;
}
