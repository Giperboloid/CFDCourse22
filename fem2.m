grid2utils;
femutils;

global grid elements u

% 1. ============================ Load grid
grid = gu_reggrid(0, 0, 1, 1, 10, 10);

Nbasis = grid.Nvert;

% elements
elements = cell(1, grid.Nelem);
for ielem = 1:grid.Nelem
	vert = cell2mat(grid.vert(grid.elem_vert{ielem})');
	elements{ielem} = D2Lin_Rect(vert, grid.elem_vert{ielem});
endfor
fem_elements_report(elements);

% 2. ============================= Analytical Functions
% Known exact solution (analytical function)
function ret = uexact(p)
	x = p(1);
	ret = sin(9*(x+0.2).^2);
endfunction

% K(x) (known analytical function)
function ret = kfun(p)
	x = p(1);
	ret = 1.0/(18*18*(x+0.2));
endfunction

% F(x) right side function (known analytical function)
function ret = ffun(p)
	x = p(1);
	ret = (x + 1.2).*uexact(x);
endfunction

% numerical solution: defined through solution vector u and basis functions phi_i
function ret = unumer(p)
	global u grid elements
	ielem = gu_find_element(grid, p);
	e = elements{ielem};
	xi = e.to_xi(p);
	ret = 0;
	for i=1:e.nbasis
		ret += u(e.ibasis(i))*e.phi{i}(xi);
	endfor
endfunction

% approximate analytic functions
fvec = [];
kvec = [];
for i=1:grid.Nelem
	pnt = grid.elem_center{i};
	kvec(i, 1) = kfun(pnt);
endfor
for i=1:grid.Nvert;
	pnt = grid.vert{i};
	fvec(i, 1) = ffun(pnt);
endfor
if Nbasis > grid.Nvert
	error("Define rhs approximation for the supplementary basis");
endif

% 3 ================================ SLE Assembly and Solution
% left hand side matrices
Mass = zeros(Nbasis, Nbasis);
Stiff = zeros(Nbasis, Nbasis);
LumpedMass = zeros(Nbasis, 1);

% element-wise assembly
for ielem = 1:grid.Nelem
	e = elements{ielem};
	for i = 1:e.nbasis
		for j = 1:e.nbasis
			% add local matrices to global mass and stiff matrices
			Mass(e.ibasis(i), e.ibasis(j)) += e.mass(i, j);
			Stiff(e.ibasis(i), e.ibasis(j)) += kvec(ielem)*e.stiff(i, j);
		endfor
		LumpedMass(e.ibasis(i)) += e.lumped_mass(i);
	endfor
endfor

% final matrices
M = Mass + Stiff;
rhs = Mass*fvec;

% Dirichlet boundary conditions
for i=1:grid.Nvert
	x = grid.vert{i}(1);
	M(i, :) = 0;
	M(i, i) = 1;
	rhs(i) = uexact(grid.vert{i});
endfor

% obtain solution
u = M\rhs;

% 4. ============================== Visualization and output
% visualization on equidistant fine grid with Nvis nodes
Nvis = 100;
x = linspace(0, 1, Nvis);
y = 0.5;
clear y_exact y_numer
for i=1:Nvis
	% get exact solution for each x entry
	y_exact(i) = uexact([x(i), y]);
	% get numerical solution for each x entry
	y_numer(i) = unumer([x(i), y]);
endfor
% plot y_exact(x), y_numer(x)
p = plot(x, y_exact, x, y_numer);
lg = legend("EXACT", "NUMER");
set(lg, "fontsize", 14);
set(p, "linewidth", 3);
set(p, "markersize", 15);

% calculate and print norms: maximum and standard deviations
Nmax = max(abs(y_exact - y_numer))
N2 = std(y_exact - y_numer)

uexact_appr = zeros(grid.Nvert, 1);
for i=1:grid.Nvert
	pnt = grid.vert{i};
	uexact_appr(i) = uexact(pnt);
endfor

NmaxA = max(abs(uexact_appr - u(1:grid.Nvert)))
N2A = sqrt(1/sum(grid.elem_volume) * dot(LumpedMass, (uexact_appr - u(1:grid.Nvert)).^2))

% vtk output
gu_save_vtk_vertdata(grid, u(1:grid.Nvert), "result.vtk")
