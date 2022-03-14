% XY-Aligned rectangle element with bilinear basis:
% 
% 4---------3
% |         |
% |         |
% 1 ------- 2
%
%  Parametric square: xi = [0, 1], eta = [0, 1]
classdef D2Lin_Rect
properties
	% number of basis functions
	nbasis = 4
	% local basis functions
	phi = { @(xi) (1 - xi(1))*(1 - xi(2));
		@(xi) xi(1)*(1-xi(2));
		@(xi) xi(1)*xi(2);
		@(xi) (1 - xi(1))*xi(2); }
	% global indices of local vertices. [1, 4] array of indices
	ibasis
	% local mass matrix. [4, 4] array
	mass
	% local stiffness matrix. [4, 4] array
	stiff
	% local lumped mass matrix. [1, 4] array
	lumped_mass
	% local transport matrix. [4, 4] array for two dimensions = {[4, 4], [4, 4]}
	tran

	_p0
	_l
end
methods
	% Constructor
	% vertices - array of vertex coordinates. Array [4, 2] (4 points with x and y coordinates)
	% ibasis - global indices of local vertices. Array [1, 4] (global indices of vertices in counter-clockwise direction)
	function ret = D2Lin_Rect(vertices, ibasis)
		assert(size(ibasis) == [1, 4]);
		assert(size(vertices) == [4, 2]);

		Lx = vertices(2, 1) - vertices(1, 1);
		Ly = vertices(4, 2) - vertices(1, 2);

		ret._p0 = vertices(1, :);
		ret._l = [Lx, Ly];

		ret.ibasis = ibasis;

		ret.mass = Lx*Ly/36 * [...
			4, 2, 1, 2;
			2, 4, 2, 1;
			1, 2, 4, 2;
			2, 1, 2, 4];

		ret.stiff = 1/(6*Lx*Ly)*[...
			 2*(Ly^2+Lx^2),  Lx^2-2*Ly^2,   -Ly^2-Lx^2,      Ly^2-2*Lx^2;
			 Lx^2-2*Ly^2,    2*(Ly^2+Lx^2),  Ly^2-2*Lx^2,   -Ly^2-Lx^2;
			-Ly^2-Lx^2,      Ly^2-2*Lx^2,    2*(Ly^2+Lx^2),  Lx^2-2*Ly^2;
			 Ly^2-2*Lx^2,   -Ly^2-Lx^2,      Lx^2-2*Ly^2,    2*(Ly^2+Lx^2)];

		ret.lumped_mass = Lx*Ly*[...
			1/4;
			1/4;
			1/4;
			1/4;
		];

		ret.tran = {...
			Ly/12*[-2, 2, 1, -1;
			       -2, 2, 1, -1;
			       -1, 1, 2, -2;
			       -1, 1, 2, -2],...
			Lx/12*[-2, -1, 1, 2;
			       -1, -2, 2, 1;
			       -1, -2, 2, 1;
			       -2, -1, 1, 2]
		}
	)

	)
	endfunction

	% convert point in physical space to point in parametric space
	%   x - Point in physical space: array[1, 2]
	%   returns xi - Point in parametric space: array[1, 2]
	function ret = to_xi(self, x)
		ret = (x-self._p0)./self._l;
	endfunction
end
end
