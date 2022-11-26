#include "test_solver.hpp"
#include <iostream>

TestSolver::TestSolver(double value){
    this->value = value;
}

void TestSolver::set_matrix(const CsrStencil& mat, const std::vector<double>& mat_values){
    
}

void TestSolver::solve(const std::vector<double>& rhs, std::vector<double>& ret) const{
    ret.resize(rhs.size());
    for (int i=0; i<rhs.size(); i++){
        ret[i]=this->value;
    }
}