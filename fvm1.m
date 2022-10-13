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
% number of faces
Nface = Nvert;
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
	ret = sin(9*(x+0.2).^2);
endfunction

% K(x) (known analytical function)
function ret = kfun(x)
	ret = 1.0/(18*18*(x+0.2));
endfunction

% F(x) right side function (known analytical function)
function ret = ffun(x)
	ret = (x + 1.2).*uexact(x);
endfunction

function ret = find_elem(x)
	global Nelem grid
	% TODO Optimize
	for ielem=1:Nelem
		if x >= grid(ielem) && x <= grid(ielem+1)
			ret = ielem;
			return
		endif
	endfor
	error("%f is out of bounds\n", x);
endfunction

% numerical solution: defined through solution vector u and basis functions phi_i
function ret = unumer(x)
	global u
	ielem = find_elem(x);
	ret = u(ielem);
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

% approximate analytic functions
fvec = [];
kvec = [];
for i=1:Nelem
	fvec(i) = ffun(center(i));
	kvec(i) = kfun(center(i));
endfor



% 3. ================================ SLE Assembly and Solution
% left hand side matrix
M = zeros(Nelem, Nelem);
% right hand side vector
rhs = zeros(Nelem, 1);

% assembly: looping through each element (or each row)
for irow = 1:Nelem
	ielem = irow;
	% --- integral 2 (mass)
	len = grid(ielem+1) - grid(ielem);
	M(irow, irow) += len;

	% --- integral 3 (right side)
	rhs(irow) += len*ffun(center(ielem));
endfor

% internal
for iface=2:Nface-1
	ielem1 = iface - 1;
	ielem2 = iface;
	irow1 = ielem1;
	irow2 = ielem2;

	k1 = kvec(ielem1);
	k2 = kvec(ielem2);
	h1 = grid(ielem1+1) - center(ielem1);
	h2 = center(ielem2) - grid(ielem2);
	kc = k1*k2*(h1 + h2)/(k1*h2 + k2*h1);

	m = kc/(center(ielem2) - center(ielem1));

	M(irow1, irow1) += m;
	M(irow2, irow2) += m;
	M(irow1, irow2) += -m;
	M(irow2, irow1) += -m;
endfor

% left side Dirichlet boundary condition
kc = kvec(1);
m = kc/(center(1) - grid(1));
M(1, 1) += m;
rhs(1) += uexact(grid(1)) * m;

% right side Dirichlet boundary condition
kc = kvec(Nelem);
m = kc/(grid(Nvert) - center(Nelem));
M(Nelem, Nelem) += m;
rhs(Nelem) += uexact(grid(Nvert)) * m;

% obtain solution
u = M\rhs;

% 4. ============================== Visualization and output
% visualization on equidistant fine grid with Nvis nodes
Nvis = 100;
x = linspace(A, B, Nvis);
clear y_exact y_numer
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

uexact_appr = zeros(Nelem, 1);
for i=1:Nelem
	uexact_appr(i) = uexact(center(i));
endfor
NmaxA = max(abs(uexact_appr - u))
N2A = sqrt(1/Nelem*sum((uexact_appr - u).^2))
