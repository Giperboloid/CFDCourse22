1;
% ===== FINITE VOLUME METHOD 2D
% ===== -d( K(x, y) * du(x, y)/dx )/dx
% ===== -d( K(x, y) * du(x, y)/dy )/dy + u(x, y) = F(x, y) ====
grid2utils;

% use some global variables for the sake of convenience
global A B Nelem Nvert grid u

%grid = gu_reggrid(0, 0, 1, 1, 30, 30);
grid = gu_build_from_gmsh_m("grid_examples/tri_1000.m");

% 1. ===================================== Input Data
% number of finite elements
Nelem = grid.Nelem;
% number of vertices
Nvert = grid.Nvert;
% number of faces
Nface = grid.Nface;

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
	global grid u
	ielem = gu_find_element(grid, p);
	ret = u(ielem);
endfunction

% approximate analytic functions
fvec = [];
kvec = [];
for i=1:Nelem
	pnt = grid.elem_center{i};
	fvec(i) = ffun(pnt);
	kvec(i) = kfun(pnt);
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
	v = grid.elem_volume(ielem);
	M(irow, irow) += v;

	% --- integral 3 (right side)
	rhs(irow) += v*fvec(ielem);
endfor

% stiffness (internal)
for iface = grid.internal_faces
	ielem1 = grid.face_elem{iface}(1);
	ielem2 = grid.face_elem{iface}(2);
	irow1 = ielem1;
	irow2 = ielem2;

	k1 = kvec(ielem1);
	k2 = kvec(ielem2);
	h1 = grid.face_elem_distance{iface}(1);
	h2 = grid.face_elem_distance{iface}(2);
	kc = k1*k2*(h1 + h2)/(k1*h2 + k2*h1);
	h = grid.face_cross_distance(iface);
	m = kc/h*grid.face_area(iface)*grid.face_cosn(iface);

	M(irow1, irow1) += m;
	M(irow2, irow2) += m;
	M(irow1, irow2) += -m;
	M(irow2, irow1) += -m;
endfor

% stiffness (Dirichlet boundaries)
for iface=grid.boundary_faces
	x1 = grid.face_center{iface}(1);
	if grid.face_elem{iface}(1) > 0
		iloc = 1;
	else
		iloc = 2;
	endif
	ielem = grid.face_elem{iface}(iloc);

	h = grid.face_cross_distance(iface);
	kc = kvec(ielem);
	m = kc/h*grid.face_area(iface)*grid.face_cosn(iface);
	M(ielem, ielem) += m;
	pnt = grid.face_center{iface};
	rhs(ielem) += uexact(pnt) * m;
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
plot(x, y_exact, x, y_numer);
lg = legend("EXACT", "NUMER");
set(lg, "fontsize", 14);

% calculate and print norms: maximum and standard deviations
Nmax = max(abs(y_exact - y_numer))
N2 = std(y_exact - y_numer)

uexact_appr = zeros(Nelem, 1);
for i=1:Nelem
	pnt = grid.elem_center{i};
	uexact_appr(i) = uexact(pnt);
endfor

NmaxA = max(abs(uexact_appr - u))
N2A = sqrt( 1/sum(grid.elem_volume) * dot(grid.elem_volume, (uexact_appr - u).^2) )

% vtk output
gu_save_vtk_elemdata(grid, u, "result.vtk")
