% Segment with linear basis
%
%  1-------2
%
%  Parametric segment: xi = [0, 1]
classdef D1Lin
properties
	% ================================ interface
	% number of basis functions
	nbasis = 2
	% local basis functions
	phi = { @(xi) 1 - xi;
		@(xi) xi }
	% global indices of local vertices. [1, 2] array of indices
	ibasis
	% local mass matrix. [2, 2] array
	mass
	% local stiffness matrix. [2, 2] array
	stiff
	% local transport matrix. [2, 2] array for one dimension = {[2, 2]}
	tran
	% local lumped mass matrix. [1, 2] array
	lumped_mass

	_p0
	_h
end
methods
	% Constructor
	% vertices - array of vertex coordinates. Array [2, 1]
	% ibasis - global indices of local vertices. Array [1, 2]
	function ret = D1Lin(vertices, ibasis)
		assert(size(ibasis) == [1, 2]);
		assert(size(vertices) == [2, 1]);

		L = vertices(2) - vertices(1);

		ret.ibasis = ibasis;
		ret._h = L;
		ret.mass = L * [1/3, 1/6; 1/6, 1/3];
		ret.stiff = 1.0/L * [ 1, -1; -1,  1];
		ret.lumped_mass = L*[1/2; 1/2];
		ret.tran = {[-1/2, 1/2; -1/2, 1/2]};

		ret._p0 = vertices(1);
	endfunction

	% convert point in physical space to point in parametric space
	function ret = to_xi(self, x)
		ret = (x - self._p0)/self._h;
	endfunction
end
end
