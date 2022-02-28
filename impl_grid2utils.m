1;
function ret = _gu_distance(a, b)
	ret = norm(a - b);
endfunction

function ret = _gu_vec_rotate(v, angle)
	c = cos(angle);
	s = sin(angle);
	ret = [c*v(1) - s*v(2), s*v(1) + c*v(2)];
endfunction

function ret = _gu_vec_cos(v1, v2)
	ret = dot(v1, v2)/norm(v1)/norm(v2);
endfunction

function ret = _gu_cross_product(a, b)
	ret = a(1)*b(2) - a(2)*b(1);
endfunction

function ret = _gu_tri_area(p1, p2, p3)
	ret = _gu_cross_product(p2 - p1, p3 - p1) / 2;
endfunction

function ret = _is_within_triangle(p1, p2, p3, p_target)
	d1 = _gu_cross_product(p1-p_target, p2-p_target);
	d2 = _gu_cross_product(p2-p_target, p3-p_target);
	d3 = _gu_cross_product(p3-p_target, p1-p_target);

	has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
	has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

	ret = !(has_neg && has_pos);
endfunction

function ielem = _gu_find_element(grid, p)
	% loop through each element
	% TODO: ------ Should be optimized
	for ie=1:grid.Nelem
		e = grid.elem_vert{ie};
		ne = grid.elem_nvert(ie);
		for i=2:ne-1
			if _is_within_triangle(grid.vert{e(1)}, grid.vert{e(i)}, grid.vert{e(i+1)}, p)
				ielem = ie;
				return
			endif
		endfor
	endfor
	error("(%f, %f) is out of bounds\n", p(1), p(2));
endfunction

function [v, cnt] = _gu_volume_center(grid)
	v = zeros(1, grid.Nelem);
	cnt = cell(1, grid.Nelem);
	for ie=1:grid.Nelem
		e = grid.elem_vert{ie};
		ne = grid.elem_nvert(ie);
		a = [];
		p0 = zeros(1, 2);
		sum_a = 0;
		for i=2:ne-1
			tri = _gu_tri_area(grid.vert{e(1)}, grid.vert{e(i)}, grid.vert{e(i+1)});
			sum_a += tri;
			p0 = p0 + tri * (grid.vert{e(1)} + grid.vert{e(i)} + grid.vert{e(i+1)})/3;
		endfor
		v(ie) = sum_a;
		cnt(ie) = p0/sum_a;
	endfor
endfunction

function vert_elem = _gu_vert_elem(grid)
	vert_elem = cell(1, grid.Nvert);

	for ie=1:grid.Nelem
		for iv=grid.elem_vert{ie}
			vert_elem{iv}(end+1) = ie;
		endfor
	endfor
endfunction

function vert_vert = _gu_vert_vert(grid)
	vert_vert = cell(1, grid.Nvert);

	for ie=1:grid.Nelem
		for iloc=1:length(grid.elem_vert{ie})
			if iloc == length(grid.elem_vert{ie})
				iloc2 = 1;
			else
				iloc2 = iloc+1;
			endif
			v1 = grid.elem_vert{ie}(iloc);
			v2 = grid.elem_vert{ie}(iloc2);
			vert_vert{v1}(end+1) = v2;
			vert_vert{v2}(end+1) = v1;
		endfor
	endfor

	for iv=1:grid.Nvert
		vert_vert{iv} = unique(vert_vert{iv});
	endfor
endfunction

function face_vert = _gu_face_vert(grid)
	face_vert = {};
	for iv1=1:grid.Nvert
		for iv2=grid.vert_vert{iv1}
			if iv1 < iv2
				face_vert{end+1} = [iv1, iv2];
			endif
		endfor
	endfor
endfunction

function vert_face = _gu_vert_face(grid)
	vert_face = cell(1, grid.Nvert);
	for iface=1:grid.Nface
		v1 = grid.face_vert{iface}(1);
		v2 = grid.face_vert{iface}(2);
		vert_face{v1}(end+1) = iface;
		vert_face{v2}(end+1) = iface;
	endfor
endfunction

function [elem_face, face_elem] = _gu_elem_face(grid)
	elem_face = cell(1, grid.Nelem);
	face_elem = cell(1, grid.Nface);
	for iface=1:grid.Nface
		face_elem{iface} = [0, 0];
	endfor

	for ie=1:grid.Nelem
		vv = grid.elem_vert{ie};
		for iloc1=1:size(vv, 2)
			v1 = vv(iloc1);
			if iloc1 == size(vv, 2)
				iloc2 = 1;
			else
				iloc2 = iloc1 + 1;
			endif
			v2 = vv(iloc2);
			if v1 > v2
				t = v2;
				v2 = v1;
				v1 = t;
				pos = false;
			else
				pos = true;
			endif
			
			for iface=1:size(grid.vert_face{v1}, 2)
				fc=grid.vert_face{v1}(iface);
				if grid.face_vert{fc}(2) == v2
					elem_face{ie}(end+1) = fc;
					if pos
						face_elem{fc}(1) = ie;
					else
						face_elem{fc}(2) = ie;
					endif
					break;
				endif
			endfor
		endfor
	endfor
endfunction

function face_center = _gu_face_center(grid)
	face_center = cell(1, grid.Nface);
	for iface=1:grid.Nface
		v1 = grid.face_vert{iface}(1);
		v2 = grid.face_vert{iface}(2);
		face_center{iface} = (grid.vert{v1} + grid.vert{v2})/2;
	endfor
endfunction

function face_area = _gu_face_area(grid);
	face_length = zeros(1, grid.Nface);
	for iface=1:grid.Nface
		v1 = grid.face_vert{iface}(1);
		v2 = grid.face_vert{iface}(2);
		face_area(iface) = _gu_distance(grid.vert{v1}, grid.vert{v2});
	endfor
endfunction

function [bnd, in] = _gu_face_types(grid)
	bnd = [];
	in = [];
	for iface=1:grid.Nface
		e = grid.face_elem{iface};
		if e(1) == 0 || e(2) == 0
			bnd(end+1) = iface;
		else
			in(end+1) = iface;
		endif
	endfor
endfunction

function _gu_save_vtk(grid, fn)
	fid = fopen(fn, 'w');
	fprintf(fid, "# vtk DataFile Version 3.0\n");
	fprintf(fid, "Grid2D\n");
	fprintf(fid, "ASCII\n");
	fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");
	% POINTS
	fprintf(fid, "POINTS %d double\n", grid.Nvert);
	for iv=1:grid.Nvert
		fprintf(fid, "%f %f 0\n", grid.vert{iv});
	endfor
	% CELLS
	fprintf(fid, "CELLS %d %d\n", grid.Nelem, grid.Nelem + sum(grid.elem_nvert));

	for ie=1:grid.Nelem
		vv = grid.elem_vert{ie};
		fprintf(fid, "%d ", size(vv, 2));
		for iv=vv
			fprintf(fid, "%d ", iv-1);
		endfor
		fprintf(fid, "\n");
	endfor

	% CELLTYPES
	types = 7*ones(grid.Nelem, 1);
	fprintf(fid, "CELL_TYPES  %d\n", grid.Nelem);
	fprintf(fid, "%d\n", types);
	
	fclose(fid);
endfunction

function _gu_vtk_add_elemdata(vec, fn)
	fid = fopen(fn, 'a');
	fprintf(fid, "CELL_DATA %d\n", length(vec));
	fprintf(fid, "SCALARS data double 1\n");
	fprintf(fid, "LOOKUP_TABLE default\n");
	fprintf(fid, "%d\n", vec);
	fclose(fid);
endfunction

function _gu_vtk_add_vertdata(vec, fn)
	fid = fopen(fn, 'a');
	fprintf(fid, "POINT_DATA %d\n", length(vec));
	fprintf(fid, "SCALARS data double 1\n");
	fprintf(fid, "LOOKUP_TABLE default\n");
	fprintf(fid, "%d\n", vec);
	fclose(fid);
endfunction

function dist = _gu_face_elem_dist(grid)
	dist = cell(1, grid.Nface);
	for iface=1:grid.Nface
		el = grid.face_elem{iface};
		d = [-1, -1];
		if el(1) > 0
			d(1) = _gu_distance(grid.face_center{iface}, grid.elem_center{el(1)});
		endif
		if el(2) > 0
			d(2) = _gu_distance(grid.face_center{iface}, grid.elem_center{el(2)});
		endif
		dist{iface} = d;
	endfor
endfunction

function [cosn, cross_distance] = _gu_face_cosn_cross_distance(grid)
	cosn = zeros(1, grid.Nface);
	cross_distance = zeros(1, grid.Nface);
	for iface=grid.internal_faces
		el = grid.face_elem{iface};
		p1 = grid.elem_center{el(1)};
		p2 = grid.elem_center{el(2)};
		cross_distance(iface) = _gu_distance(p1, p2);

		ivert1 = grid.face_vert{iface}(1);
		ivert2 = grid.face_vert{iface}(2);

		vec1 = p2 - p1;
		vec2 = grid.vert{ivert2} - grid.vert{ivert1};
		vec2 = _gu_vec_rotate(vec2, -pi/2);
		cosn(iface) = _gu_vec_cos(vec1, vec2);
	endfor

	for iface=grid.boundary_faces
		if grid.face_elem{iface}(1) > 0
			iloc = 1;
		else
			iloc = 2;
		endif
		el = grid.face_elem{iface}(iloc);

		p1 = grid.elem_center{el};
		p2 = grid.face_center{iface};
		cross_distance(iface) = _gu_distance(p1, p2);

		vec1 = p2 - p1;
		if iloc == 2
			vec1 *= -1;
		endif

		ivert1 = grid.face_vert{iface}(1);
		ivert2 = grid.face_vert{iface}(2);
		vec2 = grid.vert{ivert2} - grid.vert{ivert1};
		vec2 = _gu_vec_rotate(vec2, -pi/2);

		cosn(iface) = _gu_vec_cos(vec1, vec2);
	endfor
endfunction

function grid = _gu_assemble_grid(vert, elem_vert)
	clear grid;
	grid.Nelem = 0;
	grid.Nvert = 0;
	grid.Nface = 0;

	grid.vert = vert;
	grid.elem_vert = elem_vert;

	grid.Nelem = size(grid.elem_vert, 2);
	grid.Nvert = size(grid.vert, 2);

	grid.elem_nvert = zeros(1, grid.Nelem);
	for i=1:grid.Nelem
		grid.elem_nvert(i) = size(grid.elem_vert{i}, 2);
	endfor

	[grid.elem_volume, grid.elem_center] = _gu_volume_center(grid);
	grid.vert_elem = _gu_vert_elem(grid);
	grid.vert_vert = _gu_vert_vert(grid);
	grid.face_vert = _gu_face_vert(grid);
	grid.Nface = size(grid.face_vert, 2);
	grid.vert_face = _gu_vert_face(grid);
	[grid.elem_face, grid.face_elem] = _gu_elem_face(grid);
	grid.face_center = _gu_face_center(grid);
	grid.face_area = _gu_face_area(grid);
	[grid.boundary_faces, grid.internal_faces] = _gu_face_types(grid);
	grid.face_elem_distance = _gu_face_elem_dist(grid);
	[grid.face_cosn, grid.face_cross_distance] = _gu_face_cosn_cross_distance(grid);
	printf("== Created grid: %d elements, %d vertices\n", grid.Nelem, grid.Nvert);
endfunction
