function RR = RRes(A0, A1, V, LAMBDA)
	n = size(A0, 1);
	RR = zeros(2*n, 1);
	for i=1:2*n
		RR(i) = res(A0, A1, V(:, i), LAMBDA(i));
	endfor
	
	function val = res(A0, A1, x, l)
		num = norm((l^2)*(A1.')*x + l*A0*x + A1*x, 2);
		den = (abs(l)^2)*norm(A1, 'fro');
		den = den + abs(l)*norm(A1, 'fro');
		den = den + norm(A1, 'fro');
		den = den*norm(x, 2);
		val = num / den;
	endfunction
endfunction