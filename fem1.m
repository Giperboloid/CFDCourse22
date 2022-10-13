% ===== LINEAR FINITE ELEMENT METHOD 1D
% ===== -d( K(x) * du(x)/dx )/dx + u(x) = F(x) ====

femutils;

% use some global variables for the sake of convenience
global A B Nelem Nvert grid u center elements

% 1. ===================================== Input Data
% left boundary
A = 0;
% right boundary
B = 1;
% number of finite elements
Nelem = 10;
% number of vertices
Nvert = Nelem + 1;
% number of basis functions
Nbasis = Nvert;
% vertices of computational grid: equidistant meshing
grid = linspace(A, B, Nvert);
% solution vector
u = zeros(Nbasis, 1);
% center coordinates of finite elements
center = zeros(Nelem, 1);
for i = 1:Nelem
	center(i) = (grid(i+1) + grid(i))/2;
endfor
% elements
elements = cell(1, Nelem);
for ielem = 1:Nelem
	vert = [grid(ielem); grid(ielem+1)];
	elements{ielem} = D1Lin(vert, [ielem, ielem+1]);
endfor
fem_elements_report(elements);


% 2. ============================= Analytical Functions
% Known exact solution (analytical function)
function ret = uexact(x)
	%ret = 3*x^5 -1* x^4 -1*x^3 - 3*x^2 + x;
	ret = sin(9*(x+0.2).^2);
endfunction

% K(x) (known analytical function)
function ret = kfun(x)
	%ret = 1;
	ret = 1.0/(18*18*(x+0.2));
endfunction

% F(x) right side function (known analytical function)
function ret = ffun(x)
	%ret = -(60*x^3 - 12*x^2 - 6*x - 6) + uexact(x);
	ret = (x + 1.2)*uexact(x);
endfunction

function ret = find_elem(x)
	global grid Nelem Nvert
	if x < grid(1) || x > grid(Nvert)
		error("%f is out of bounds\n", x);
	endif
	ret = min(Nelem, max(1, floor((x - grid(1))/(grid(2) - grid(1))) + 1));
endfunction

% numerical solution: defined through solution vector u and basis functions phi_i
function ret = unumer(x)
	global u elements
	ielem = find_elem(x);
	e = elements{ielem};
	xi = e.to_xi(x);
	ret = 0;
	for i=1:e.nbasis
		ret += u(e.ibasis(i))*e.phi{i}(xi);
	endfor
endfunction

% approximate analytic functions
fvec = [];
kvec = [];
for i=1:Nelem
	kvec(end+1, 1) = kfun(center(i));
endfor
for i=1:Nvert
	fvec(end+1, 1) = ffun(grid(i));
endfor
assert(size(fvec, 1) == Nbasis);

% 3.2 ================================ SLE Assembly and Solution
% left hand side matrices
Mass = zeros(Nbasis, Nbasis);
Stiff = zeros(Nbasis, Nbasis);

% element-wise assembly
for ielem = 1:Nelem
	e = elements{ielem};
	for i = 1:e.nbasis
	for j = 1:e.nbasis
		Mass(e.ibasis(i), e.ibasis(j)) += e.mass(i, j);
		Stiff(e.ibasis(i), e.ibasis(j)) += kvec(ielem)*e.stiff(i, j);
	endfor
	endfor
endfor

% final matrices
M = Mass + Stiff;
rhs = Mass*fvec;

% left side Dirichlet boundary condition
M(1, :) = 0;
M(1, 1) = 1;
rhs(1) = uexact(grid(1));

% right side Dirichlet boundary condition
M(Nvert, :) = 0;
M(Nvert, Nvert) = 1;
rhs(Nvert) = uexact(grid(Nvert));

% obtain solution
u = M\rhs;

% 4. ============================== Visualization and output
% visualization on equidistant fine grid with Nvis nodes
Nvis = 1000;
x = linspace(A, B, Nvis);
clear y_exact y_numer
for i=1:Nvis
	% get exact solution for each x entry
	y_exact(i) = uexact(x(i));
	% get numerical solution for each x entry
	y_numer(i) = unumer(x(i));
endfor
% plot y_exact(x), y_numer(x)
p = plot(x, y_exact, x, y_numer, grid, u(1:Nvert), "x");
lg = legend("EXACT", "NUMER");
set(lg, "fontsize", 14);
set(p, "linewidth", 3);
set(p, "markersize", 15);

% calculate and print norms: maximum and standard deviations
Nmax = max(abs(y_exact - y_numer))
N2 = std(y_exact - y_numer)

uexact_appr = zeros(Nvert, 1);
for i=1:Nvert
	uexact_appr(i) = uexact(grid(i));
endfor
NmaxA = max(abs(uexact_appr - u(1:Nvert)))
N2A = sqrt(1/Nelem*sum((uexact_appr - u(1:Nvert)).^2))
