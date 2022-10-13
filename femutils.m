1;
function fem_elements_report(elements)
	all_classes = {
		"Unknown", ...
		"D1Lin", ...
		"D1Sqr", ...
		"D1Hermite", ...
		"D2Lin_Quad", ...
		"D2Lin_Rect" ...
	};
	nclasses = length(all_classes);
	res = containers.Map();
	for c=1:nclasses
		res(all_classes{c}) = 0;
	endfor
	for ielem=1:size(elements, 2)
		el = elements{ielem};
		r = "Unknown";
		for c=1:nclasses
			if isa(el, all_classes{c})
				r = all_classes{c};
				break;
			endif
		endfor
		res(r) += 1;
	endfor
	printf("== Assembled %d finite elements:\n", size(elements, 2));
	for c=1:nclasses
		key = all_classes{c};
		num = res(key);
		if num > 0
			printf("    %s - %d\n", key, num);
		endif
	endfor
endfunction

function ret1 = fem_algebraic_upwinding(Tran)
	n = size(Tran, 1);
	for i = 1:n
	for j = i+1:n
		if Tran(i, j) > 0
			diff = Tran(i, j);
			Tran(i, i) += diff;
			Tran(j, j) += diff;
			Tran(i, j) -= diff;
			Tran(j, i) -= diff;
		endif
	endfor
	endfor
	ret1 = Tran;
endfunction
