% Segment with square basis
%
%  1---3---2
%
classdef D1Sqr
properties
	nbasis = 3
	phi = { @(xi) 2*xi*xi - 3*xi + 1;
		@(xi) 2*xi*xi - xi;
		@(xi) -4*xi*xi + 4*xi;}
	ibasis
	mass
	stiff
	lumped_mass
	tran

	_p0
	_h
end
methods
	function ret = D1Sqr(vertices, ibasis)
		assert(size(ibasis) == [1, 3]);
		assert(size(vertices) == [2, 1]);

		L = vertices(2) - vertices(1);

		ret.ibasis = ibasis;
		ret._p0 = vertices(1);
		ret._h = L;
		ret.mass = L/30 * [... 
			4, -1,  2;
			-1, 4,  2;
			2,  2, 16];
		ret.stiff = 1/(3*L) * [...
			7,   1, -8;
			1,   7, -8;
			-8, -8, 16];

		ret.lumped_mass = L * [1/6; 1/6; 2/3];
		ret.tran = [...
			-1/2, -1/6, 2/3;
			1/6, 1/2, -2/3;
			-2/3, 2/3, 0];
	endfunction

	function ret = to_xi(self, x)
		ret = (x - self._p0)/self._h;
	endfunction
end
end
