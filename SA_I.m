function [V_I, LAMBDA_I] = SA_I (A0, A1, verbose = false)
	verbose && fprintf('--- SA_I ---\n');
	n = size(A0, 1);
	
	verbose && fprintf('- Step 1/5'); tic;
	% Step 1: Form the pair (K, N) as in (2.4)	
	K = [A0 A1.'-A1; A1-A1.' A0];
	N = blkdiag(-A1, -A1.');
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 2/5'); tic;
	% Step 2: Reduce (K, N) to block upper triangular forms as in (2.5) using
	%	unitary transformations
	[Kn, Nn, ~, Z] = rbutf(K, N);
	K11 = Kn(1:n, 1:n);
	N11 = Nn(1:n, 1:n);
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 3/5'); tic;
	% Step 3: Compute eigenpairs of (K11, N11)
	[~, ~, ~, ~, Y, ~, Mu] = qz(K11, N11);
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 4/5'); tic;
	% Step 4: Compute zi
	J = [zeros(n) eye(n); -eye(n) zeros(n)];
	Z = (J.') * Z * [Y; zeros(n)];
	Z1 = Z(1:n, :);
	Z2 = Z(n+1:2*n, :);
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 5/5'); tic;
	% Step 5: Compute vi, 1/vi, xi1 and xi2
	Nu = complex(zeros(n,1), zeros(n,1));
	invNu = complex(zeros(n,1), zeros(n,1));
	X1 = complex(zeros(n, n), zeros(n, n));
	X2 = complex(zeros(n, n), zeros(n, n));
	for i=1:n
		r = roots([1 -Mu(i) 1]);
		Nu(i) = r(1);
		invNu(i) = r(2);
		X1(:, i) = Z1(:, i) + invNu(i) * Z2(:, i);
		X2(:, i) = Z1(:, i) + Nu(i) * Z2(:, i);
	endfor
	verbose && fprintf(': %f s\n', toc());
	
	V_I = [X1 X2];
	LAMBDA_I = [Nu; invNu];
endfunction