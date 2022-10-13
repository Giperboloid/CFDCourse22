1;
% ===== LINEAR FINITE ELEMENT METHOD 1D
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
u = zeros(Nvert, 1);
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
	ret = (x + 1.2)*uexact(x);
endfunction

function ret = find_elem(x)
	global grid Nelem
	if x < grid(1) || x > grid(end)
		error("%f is out of bounds\n", x);
	endif
	ret = min(Nelem, max(1, floor((x - grid(1))/(grid(2) - grid(1))) + 1));
endfunction

% numerical solution: defined through solution vector u and basis functions phi_i
function ret = unumer(x)
	global u
	ielem = find_elem(x);
	ret = u(ielem)*phi(ielem, x) + u(ielem+1)*phi(ielem+1, x);
endfunction

% basis function for i-th element
function ret = phi(i, x)
	global grid Nvert
	if i > 1 && x >= grid(i-1) && x <= grid(i)
		ret = (x - grid(i-1))/(grid(i) - grid(i-1));
	elseif i < Nvert && x >= grid(i) && x <= grid(i+1)
		ret = 1 - (x - grid(i))/(grid(i+1) - grid(i));
	else
		ret = 0;
	endif
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

% 3. ================================ SLE Assembly and Solution
% left hand side matrices
Mass = zeros(Nvert, Nvert);
Stiff = zeros(Nvert, Nvert);

%elements = cell(Nelem, 1);
%for ielem = 1:Nelem
%        ivert1 = ielem;
%        ivert2 = ielem+1;
%        h = grid(ivert2) - grid(ivert1);
%        k = kvec(ielem);

%        elements{ielem}.nbasis = 2;
%        elements{ielem}.ibasis = [ielem, ielem+1];
%        elements{ielem}.mass = [h/3, h/6;h/6, h/3];
%        elements{ielem}.stiff = [k/h, -k/h; -k/h, k/h];
%endfor

% element-wise assembly
for ielem = 1:Nelem
	%elem = elements{ielem};
	%for i = 1:elem.nbasis
	%for j = 1:elem.nbasis
	%        gi = elem.ibasis(i);
	%        gj = elem.ibasis(j);
	%        Mass(gi, gj) += elem.mass(i, j);
	%        Stiff(gi, gj) += elem.stiff(i, j);
	%endfor
	%endfor
	ivert1 = ielem;
	ivert2 = ielem+1;
	h = grid(ivert2) - grid(ivert1);
	k = kvec(ielem);

	% Mass
	Mass(ivert1, ivert1) += h/3;
	Mass(ivert2, ivert2) += h/3;
	Mass(ivert1, ivert2) += h/6;
	Mass(ivert2, ivert1) += h/6;

	% Stiff
	Stiff(ivert1, ivert1) += k/h;
	Stiff(ivert2, ivert2) += k/h;
	Stiff(ivert1, ivert2) -= k/h;
	Stiff(ivert2, ivert1) -= k/h;
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
p = plot(x, y_exact, x, y_numer);
lg = legend("EXACT", "NUMER");
set(lg, "fontsize", 14);
set(p, "linewidth", 3);

% calculate and print norms: maximum and standard deviations
Nmax = max(abs(y_exact - y_numer))
N2 = std(y_exact - y_numer)

uexact_appr = zeros(Nvert, 1);
for i=1:Nvert
	uexact_appr(i) = uexact(grid(i));
endfor
NmaxA = max(abs(uexact_appr - u))
N2A = sqrt(1/Nelem*sum((uexact_appr - u).^2))
