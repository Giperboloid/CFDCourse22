1;
% ===== FINITE VOLUME METHOD 1D
% ===== -d( K(x) * du(x)/dx )/dx + u(x) = F(x) ====

% use some global variables for the sake of convenience
global A B Nelem Nvert grid u center

% 1. ===================================== Input Data
% left boundary
A = 0;
% right boundary
B = 1;
% number of finite elements
Nelem = 10;
% number of vertices
Nvert = Nelem + 1;
% vertices of computational grid: equidistant meshing
grid = linspace(A, B, Nvert);
% solution vector
u = zeros(Nelem, 1);
% center coordinates of finite elements
center = zeros(Nelem, 1);
for i = 1:Nelem
	center(i) = (grid(i+1) + grid(i))/2;
endfor

% 2. ============================= Analytical Functions
% Known exact solution (analytical function)
function ret = uexact(x)
	ret = x+1;
endfunction

% K(x) (known analytical function)
function ret = kfun(x)
	ret = 1;
endfunction

% F(x) right side function (known analytical function)
function ret = ffun(x)
	ret = uexact(x);
endfunction

% numerical solution: defined through solution vector u and basis functions phi_i
function ret = unumer(x)
	global u Nelem
	ret = 0;
	for i=1:Nelem
		ret += u(i)*phi(i, x);
	endfor
endfunction

% basis function for i-th element
function ret = phi(i, x)
	global grid
	if x >= grid(i) && x <= grid(i+1)
		ret = 1;
	else
		ret = 0;
	endif
endfunction


% 3. ================================ SLE Assembly and Solution
% left hand side matrix
M = zeros(Nelem, Nelem);
% right hand side vector
rhs = zeros(Nelem, 1);

% assembly: looping through each element (or each row)
for irow = 1:Nelem
	ielem = irow;
	% --- integral 1 (stiffness)
	% left neighbour
	if ielem > 1
		% internal
		m = kfun(grid(ielem))/(center(ielem) - center(ielem-1));
		M(irow, irow) += m;
		M(irow, irow-1) += -m;
	else
		% left side Dirichlet boundary condition
		m = kfun(grid(1))/(center(1) - grid(1));
		M(irow, irow) += m;
		rhs(irow) += uexact(grid(1)) * m;
	endif
	% right neighbour
	if ielem < Nelem
		% internal
		m = kfun(grid(ielem+1))/(center(ielem+1) - center(ielem));
		M(irow, irow) += m;
		M(irow, irow+1) += -m;
	else
		% right side Dirichlet boundary condition
		m = kfun(grid(Nvert))/(grid(Nvert) - center(Nelem));
		M(irow, irow) += m;
		rhs(irow) += uexact(grid(Nvert)) * m;
	endif

	% --- integral 2 (mass)
	len = grid(ielem+1) - grid(ielem);
	M(irow, irow) += len;

	% --- integral 3 (right side)
	rhs(irow) += len*ffun(center(ielem));
endfor

% obtain solution
u = M\rhs;

% 4. ============================== Visualization and output
% visualization on equidistant fine grid with Nvis nodes
Nvis = 1000;
x = linspace(A, B, Nvis);
for i=1:Nvis
	% get exact solution for each x entry
	y_exact(i) = uexact(x(i));
	% get numerical solution for each x entry
	y_numer(i) = unumer(x(i));
endfor
% plot y_exact(x), y_numer(x)
plot(x, y_exact, x, y_numer);
lg = legend("EXACT", "NUMER");
set(lg, "fontsize", 14);

% calculate and print norms: maximum and standard deviations
Nmax = max(abs(y_exact - y_numer))
N2 = std(y_exact - y_numer)
