#include "slae/jacobi_solver.hpp"
#include "slae/csrmat.hpp"
#include <string>

double vec_vec_residual(const std::vector<double>& v1, const std::vector<double>& v2){
	double s = 0.0;
	for (int i=0;i<v1.size(); i++){
		s += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	return sqrt(s);
}

bool test_1(double eps){
    /* {{5.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {0.0, 0.0, 6.0}}.{1.0, 2.0, 3.0} == {5.0, 5.0, 18.0}
        addres = [0, 1, 3, 4]
        columns = [1, 1, 2, 3]
        values = [5.0, 1.0, 2.0, 6.0]*/
    std::vector<std::set<int>> stencil_set{{0}, {0, 1}, {2}};
    std::vector<double> vals{5.0, 2.0, 1.0, 6.0}; 
    std::vector<double> rhs{5.0, 5.0, 18.0}; 

    std::vector<double> answer{1.0, 2.0, 3.0}; 
    std::vector<double> result; 
    
    CsrStencil st = CsrStencil::build(stencil_set);
    JacobiSolver solver(0.0001, 1000);
    solver.set_matrix(st, vals);
    solver.solve(rhs, result);
    double err = vec_vec_residual(result, answer);
    
    return err < eps;
}

bool test_2(double eps){
    /* {{3.0, 2.3, 0.0, 0.0}, {0.0, 3.5, 0.0, 0.0}, {0.0, 0.0, 8.0, 7.0}, {-1.0, 0.0, 3.0, 15.0}}.{2., 0., 9., -1.} = {6., 0., 65., 10.}
        addres = [0, 2, 3, 7, 8]
        values = [3.0, 2.3, 3.5, 8.0, 7.0, -1.0, 3.0, 15.0]
        columns = [1, 2, 2, 3, 4, 1, 3, 4]*/
    std::vector<std::set<int>> stencil_set{{0, 1}, {1}, {2, 3}, {0, 2, 3}};
    std::vector<double> vals{3.0, 2.3, 3.5, 8.0, 7.0,15.0, -1.0, 3.0}; 
    std::vector<double> rhs{6.0, 0.0, 65.0, 10.0}; 

    std::vector<double> answer{2.0, 0.0, 9.0, -1.0}; 
    std::vector<double> result;  
    
    CsrStencil st = CsrStencil::build(stencil_set);
    JacobiSolver solver(0.0001, 1000);
    solver.set_matrix(st, vals);
    solver.solve(rhs, result);
    double err = vec_vec_residual(result, answer);
    
    return err < eps;
}

bool test_3(){
    /* {{3.0, 2.3, 0.0}, {0.0, 0.1, 1.0}, {16.0, -20.0, 8.0}}.x == {0, 0, 0}
        values = [3, 2.3, 0.0, 1.0, 8.0, 16.0, -20.0]
        */
    std::vector<std::set<int>> stencil_set{{0, 1}, {2}, {0, 1, 2}};
    std::vector<double> vals{3.0, 2.3, 0.1, 1.0, 8.0, 16.0, -20.0}; 
    std::vector<double> rhs{0.0, 0.0, 0.0}; 

    std::vector<double> result;  
    
    CsrStencil st = CsrStencil::build(stencil_set);
    JacobiSolver solver(0.0001, 1000);
    solver.set_matrix(st, vals);
    solver.solve(rhs, result);
    
    return true;
}

bool test_4(){
    /* {{3.0, 2.3, 0.0}, {0.0, 0.0, 1.0}, {16.0, -20.0, 8.0}}.x == {0, 0, 0}
        values = [3, 2.3, 0.0, 1.0, 8.0, 16.0, -20.0]
        */
    std::vector<std::set<int>> stencil_set{{0, 1}, {2}, {0, 1, 2}};
    std::vector<double> vals{3.0, 2.3, 0.0, 1.0, 8.0, 16.0, -20.0}; 
    std::vector<double> rhs{0.0, 0.0, 0.0}; 

    std::vector<double> result;  
    
    CsrStencil st = CsrStencil::build(stencil_set);
    JacobiSolver solver(0.0001, 1000);
    try
    {
        solver.set_matrix(st, vals);
        return false;
    }
    catch(const std::exception& e)
    {
        return true;
    }
}

bool test_5(){
    /* {{5.0, 0.0, 0.0}, {1.0, 2.0, 0.0}, {0.0, 0.0, 6.0}}.{1.0, 2.0, 3.0} == {5.0, 5.0, 18.0}
        addres = [0, 1, 3, 4]
        columns = [1, 1, 2, 3]
        values = [5.0, 1.0, 2.0, 6.0]*/
    std::vector<std::set<int>> stencil_set{{0}, {0, 1}, {2}};
    std::vector<double> vals{5.0, 2.0, 1.0, 6.0}; 
    std::vector<double> rhs{5.0, 5.0, 18.0, 10.0}; // Неверный размер

    std::vector<double> answer{1.0, 2.0, 3.0}; 
    std::vector<double> result; 
    
    CsrStencil st = CsrStencil::build(stencil_set);
    JacobiSolver solver(0.0001, 1000);
    try
    {
        solver.set_matrix(st, vals);
        solver.solve(rhs, result);
        return false;
    }
    catch(const std::exception& e)
    {
        return true;
    }
}

void print_test_result(bool result, std::string test_name){
    if (result)
        std::cout << "Test \"" << test_name << "\" completed\n";
    else
        std::cout << "Test \"" << test_name << "\" failed\n";
}


int main(){
    print_test_result(test_1(0.001), "Jacoby Solver: test_1");
    print_test_result(test_2(0.001), "Jacoby Solver: test_2");
    std::cout << "Warn test-------\n";
    test_3();
    std::cout << "Warn test end---\n";
    print_test_result(test_4(), "Jacoby Solver: test_4");
    print_test_result(test_5(), "Jacoby Solver: test_5");
    return 0;
}