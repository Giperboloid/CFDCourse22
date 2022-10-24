#ifndef DBG_PRINTER_HPP
#define DBG_PRINTER_HPP

#include <vector>
#include "slae/csrmat.hpp"

namespace dbg{

void print_matrow(int irow, const CsrStencil& stencil, const std::vector<double>& mat);

}

#endif
