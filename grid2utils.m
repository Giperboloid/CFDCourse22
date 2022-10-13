1;
% Functions for working with 2D grids
% Grid structure:
%        Grid{
%                -- Data sizes
%                Nvert                   - Number of vertices
%                Nface                   - Number of faces
%                Nelem                   - Number of elements
%
%                -- Geometry data
%                vert                    - vertex coordinates, [Nvert, 2] array
%                face_center             - face centers, [Nface, 2] array
%                elem_center             - elements centers, [Nelem, 2] array

%                elem_volume             - volumes of elements, [Nelem] array
%                face_area               - face areas, [Nface] array
%                face_elem_distance      - distances from face center to connected element centers, [Nface, 2] array.
%                face_cross_distance     - distance between
%                                              left and right element centers (for internal face) 
%                                              element center and face center (for boundary face)
%                                          [Nface] array
%                face_cosn               - cosine of angle between
%                                              face right normal and vector from left element center to right element center (for internal face)
%                                              face external normal and vector between connected element center and face center (for boundry face)
%                                          [Nface] array
%
%                -- Connectivity tables
%                elem_vert               - Element->Vertex, [Nelem, {var}] array. In counterclockwise direction.
%                vert_elem               - Vertex->Element, [Nvert, {var}] array
%                vert_vert               - Vertex->Vertex, [Nvert, {var}] array
%                face_vert               - Face->Vertex, [Nface, 2] array
%                vert_face               - Vertex->Face, [Nvert, {var}] array
%                elem_face               - Element->Face, [Nelem, {var}] array
%                face_elem               - Face->Element, [Nface, 2] array.  For each face has the form [left_element, right_element]. 
%                                          For boundary faces corresponding values are zero.
%                boundary_faces          - Indices of boundary faces, [{var}] array
%                internal_faces          - Indices of internal faces, [{var}] array
%
%
%        };

% Constructors:
%        gu_reggrid(x0, y0, x1, y1, nx, ny)  - construct structured grid
%                x0, y0 - bottom left corner
%                x1, y1 - top rigth corner
%                nx, ny - partitions in x and y directions
%        gu_unstructured(vert, elem_vert) - build from custom vertices table and element-vertex connectivity
%                vert - vertices table. In form [{0, 0}, {0, 1}, {0, 2}, ....]
%                elem_vert - element-vertex connectivity. In counterclockwise direction. Has the form [{1, 2, 3}, {4, 7, 5, 2}, ...]
%        gu_build_from_gmsh_m(fn) - build from *.m file exported from gmsh
%                fn - *.m file path
% 
% Procedures:
%        gu_find_element(grid, p) - finds index of the element containing point p
%        gu_save_vtk(grid, fn) - exports grid to vtk format
%        gu_save_vtk_elemdata(grid, vec, fn) exports grid with corresponding cell data array to vtk format
%        gu_save_vtk_vertdata(grid, vec, fn) exports grid with corresponding vertex data array to vtk format

impl_grid2utils
function grid = gu_reggrid(x0, y0, x1, y1, nx, ny)
	vert = cell(1, (nx+1)*(ny+1)); 
	elem_vert = cell(1, nx*ny);

	hx = (x1 - x0)/nx;
	hy = (y1 - y0)/ny;
	for j=1:ny+1
	for i=1:nx+1
		k = (j-1)*(nx+1) + i;
		vert{k} = [x0 + hx*(i-1), y0 + hy*(j-1)];
	endfor
	endfor

	for j=1:ny
	for i=1:nx
		k = (j-1)*nx + i;
		elem_vert{k} = [ ...
			(j-1)*(nx+1) + i    ...
			(j-1)*(nx+1) + i+1  ...
			(j-0)*(nx+1) + i+1  ...
			(j-0)*(nx+1) + i    ...
		];
	endfor
	endfor

	grid = _gu_assemble_grid(vert, elem_vert);
endfunction

function grid = gu_unstructured(vert, elem_vert)
	grid = _gu_assemble_grid(vert, elem_vert);
endfunction

function grid = gu_build_from_gmsh_m(fn)
	source(fn);
	vert = {};
	for v=msh.POS'
		vert{end+1} = v'(1:2);
	endfor

	elem_vert = {};
	if isfield(msh, "TRIANGLES")
		for p = msh.TRIANGLES'
			elem_vert{end+1} = [p(1), p(2), p(3)];
		endfor
	endif
	if isfield(msh, "QUADS")
		for p = msh.QUADS'
			elem_vert{end+1} = [p(1), p(2), p(3), p(4)];
		endfor
	endif

	grid = _gu_assemble_grid(vert, elem_vert);
endfunction

function ielem = gu_find_element(grid, p)
	ielem = _gu_find_element(grid, p);
endfunction

function gu_save_vtk(grid, fn)
	_gu_save_vtk(grid, fn);
endfunction

function gu_save_vtk_elemdata(grid, vec, fn)
	_gu_save_vtk(grid, fn);
	_gu_vtk_add_elemdata(vec, fn);
endfunction

function gu_save_vtk_vertdata(grid, vec, fn)
	_gu_save_vtk(grid, fn);
	_gu_vtk_add_vertdata(vec, fn);
endfunction
