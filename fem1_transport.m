% ===== LINEAR FINITE ELEMENT METHOD 1D
% ===== dU/dt + v * dU/dx = 0
% ===== Implicit solver

femutils;

% use some global variables for the sake of convenience
global A B Nelem Nvert grid u center elements v

% 1. ===================================== Input Data
% problem dimension
ndim = 1;
% left boundary
A = 0;
% right boundary
B = 1;
% number of finite elements
Nelem = 30;
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
	conn = [ielem, ielem+1];
	% patch to connectivity due to peridicity
	if ielem == Nelem
		conn(2) = 1;
	endif
	elements{ielem} = D1Lin(vert, conn);
endfor
fem_elements_report(elements);


% 2. ============================= Analytical Functions
% velocity for one dimension
v = [1];

% Exact solution at t=0
function ret = uexact0(x)
	% periodic in [-0.5, 0.5]
	while x < -0.5
		x += 1;
	endwhile
	while x > 0.5
		x -= 1;
	endwhile

	ret=exp(-(10*x)^2);
	
	%if x > -0.1 && x < 0.1
	%        ret = 1;
	%else
	%        ret = 0;
	%endif
endfunction

% Known exact solution at t
function ret = uexact(x, t)
	global v
	ret = uexact0(x-v*t);
endfunction

% F(x) right side function (known analytical function)
function ret = ffun(x)
	ret = 0;
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

% 3.2 ================================ SLE Assembly and Solution
% left hand side matrices
Mass = zeros(Nbasis, Nbasis);
Tran = zeros(Nbasis, Nbasis);

% element-wise assembly
for ielem = 1:Nelem
	e = elements{ielem};
	for i = 1:e.nbasis
	for j = 1:e.nbasis
		Mass(e.ibasis(i), e.ibasis(j)) += e.mass(i, j);
		% = vx*grad_x(u) + vy*grad_y(u) + vz*grad_z(u)
		for idim = 1:ndim
			Tran(e.ibasis(i), e.ibasis(j)) += v(idim)*e.tran{idim}(i, j);
		endfor 
	endfor
	endfor
endfor

% algebraic upwinding
Tran = fem_algebraic_upwinding(Tran);

% time dependant parameters
tend = 0.5;
tau = 0.01;

% solution at t=0
u = zeros(Nbasis, 1);
for ivert = 1:Nvert
	u(ivert) = uexact0(grid(ivert));
endfor

% final matrices
M = 1/tau * Mass + Tran;

% Supplement zero row (periodicity)
M(Nvert, Nvert) = 1;
M(Nvert, 1) = -1;

% obtain solution in time loop
t = 0;
while t < tend
	% set uold
	t += tau;
	uold = u;

	% assemble rhs
	rhs = (1/tau)*Mass*uold;
	% last row (peridicity)
	rhs(Nvert) = 0;

	% solve matrix
	printf("t = %f\n", t);
	u = M\rhs;
endwhile

% 4. ============================== Visualization and output
% visualization on equidistant fine grid with Nvis nodes
Nvis = 1000;
x = linspace(A, B, Nvis);
clear y_exact y_numer
for i=1:Nvis
	% get exact solution for each x entry
	y_exact(i) = uexact(x(i), t);
	% get numerical solution for each x entry
	y_numer(i) = unumer(x(i), t);
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
	uexact_appr(i) = uexact(grid(i), t);
endfor
NmaxA = max(abs(uexact_appr - u(1:Nvert)))
N2A = sqrt(1/Nelem*sum((uexact_appr - u(1:Nvert)).^2))
