% Arbitrary quadrangle element with bilinear basis:
%
% 4---------3
% |         |
% |         |
% 1 ------- 2
%
classdef D2Lin_Quad
properties
	nbasis = 4
	phi = { @(xi) (1 - xi(1))*(1 - xi(2));
		@(xi) xi(1)*(1-xi(2));
		@(xi) xi(1)*xi(2);
		@(xi) (1 - xi(1))*xi(2); }
	ibasis
	mass
	stiff
	lumped_mass

	_p0
	_inv_jac
end
methods
	function ret = D2Lin_Quad(vertices, ibasis)
		assert(size(ibasis) == [1, 4]);
		assert(size(vertices) == [4, 2]);

		p0 = vertices(1, :);

		v2 = vertices(2, :) - p0;
		v3 = vertices(3, :) - p0;
		v4 = vertices(4, :) - p0;
		x2 = v2(1); x3 = v3(1); x4 = v4(1);
		y2 = v2(2); y3 = v3(2); y4 = v4(2);

		jac = @(xi, eta) [x2*(1-eta) + x3*eta - x4*eta, x2*(-xi) + x3*xi + x4*(1-xi);
				  y2*(1-eta) + y3*eta - y4*eta, y2*(-xi) + y3*xi + y4*(1-xi);];
		dxi = {...
			@(xi, eta) eta - 1;
			@(xi, eta) 1 - eta;
			@(xi, eta) eta;
			@(xi, eta) -eta;
		};
		deta = {...
			@(xi, eta) xi - 1;
			@(xi, eta) -xi;
			@(xi, eta) xi;
			@(xi, eta) 1 - xi;
		};


		px = [0.774596669241483e0, -0.774596669241483e0, 0.774596669241483e0, ...
		      -0.774596669241483e0, 0.774596669241483e0, -0.774596669241483e0, ...
		      0.0e0, 0.0e0, 0.0e0];

		py = [0.774596669241483e0, 0.774596669241483e0, -0.774596669241483e0, ...
		      -0.774596669241483e0, 0.0e0, 0.0e0, ...
		      0.774596669241483e0, -0.774596669241483e0, 0.0e0];

		w = [0.308641975308642e0, 0.308641975308642e0, 0.308641975308642e0, ...
		     0.308641975308642e0, 0.493827160493827e0, 0.493827160493827e0, ...
		     0.493827160493827e0, 0.493827160493827e0, 0.790123456790123e0];

		for i=1:4
		for j=i:4
			mass_val = zeros(9, 1);
			stiff_val = zeros(9, 1);
			for k=1:9
				xi = (px(k)+1)/2;
				eta = (py(k)+1)/2;

				J = jac(xi, eta);
				modJ = det(J);
				d_i_dxi  = dxi{i}(xi, eta);
				d_j_dxi  = dxi{j}(xi, eta);
				d_i_deta = deta{i}(xi, eta);
				d_j_deta = deta{j}(xi, eta);

				mass_val(k) = ret.phi{i}([xi, eta])*ret.phi{j}([xi, eta])*modJ;

				dx_i = J(2, 2)*d_i_dxi - J(2, 1)*d_i_deta;
				dx_j = J(2, 2)*d_j_dxi - J(2, 1)*d_j_deta;
				dy_i = -J(1, 2)*d_i_dxi + J(1, 1)*d_i_deta;
				dy_j = -J(1, 2)*d_j_dxi + J(1, 1)*d_j_deta;
				stiff_val(k) = (dx_i*dx_j + dy_i*dy_j)/modJ;
			endfor
			ret.mass(i, j) = ret.mass(j, i) = dot(mass_val, w)/4;
			ret.stiff(i, j) = ret.stiff(j, i) = dot(stiff_val, w)/4;
		endfor
		endfor

		ret.lumped_mass = [
			sum(ret.mass(1, :));
			sum(ret.mass(2, :));
			sum(ret.mass(3, :));
			sum(ret.mass(4, :));
		];

		ret.ibasis = ibasis;
		ret._p0 = p0;
		ret._inv_jac = @(xi) inv(jac(xi(1), xi(2)));
	endfunction

	function ret = to_xi(self, x)
		x -= self._p0;
		ret = (self._inv_jac(x) * x')';
	endfunction
end
end
