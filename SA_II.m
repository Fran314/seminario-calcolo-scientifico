function [V_II, LAMBDA_II] = SA_II (A0, A1, verbose = false)
	verbose && fprintf('--- SA_II ---\n');
	n = size(A0, 1);
	
	verbose && fprintf('- Step 1/5'); tic;
	% Step 1: Form the pair (Kt, Nt) as in (3.3)
	Z = [A1 A1; (A0-A1.') A1];
	Kt = [zeros(2*n) (Z.'-Z); (Z-Z.') zeros(2*n)];
	Nt = [-Z zeros(2*n); zeros(2*n) (-Z.')];
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 2/5'); tic;
	% Step 2: Reduce (Kt, Nt) to block upper triangular forms as in (3.7) using
	%	unitary transformations of (3.4)-(3.6)
	[Ktn, Ntn, ~, Va] = rbutf(Kt, Nt);
	
	P2n = complex(zeros(2*n));
	zr = complex(zeros(2*n));
	for i=1:n
		P2n(i, 2*i - 1) = 1 +0i;
		P2n(n+i, 2*i) = 1 + 0i;
	endfor
	V = Va * [P2n.' zr; zr P2n.'];
	Ktn = [P2n zr; zr P2n] * Ktn * [P2n.' zr; zr P2n.'];
	Ntn = [P2n zr; zr P2n] * Ntn * [P2n.' zr; zr P2n.'];
	R1 = Ktn(1:n, n+1:2*n);
	R2 = Ktn(n+1:2*n, 1:n);
	R3 = Ntn(1:n, 1:n);
	R4 = Ntn(n+1:2*n, n+1:2*n);
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 3/5'); tic;
	% Step 3: Compute eigenpairs {(γ i , y i )} ni=1 of (R1 * R4−1 * R2, R3 ) 
	%	in (3.8) by the periodic QZ algorithm
	[~, ~, ~, ~, Y, ~, Gamma] = qz(R1 * inv(R4) * R2, R3);
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 4/5'); tic;
	% Step 4: Compute Zt
	Mu = sqrt(Gamma);
	Eta = (inv(R1) * R3 * Y) .* (Mu.');
	Zt = V * [Y; Eta; complex(zeros(n)); complex(zeros(n))];
	Zt1 = Zt(1:2*n, :);
	Zt2 = Zt(2*n+1:4*n, :);
	verbose && fprintf(': %f s\n', toc());
	
	verbose && fprintf('- Step 5/5'); tic;
	% Step 5: Compute vi, 1/vi, xi1 and xi2
	Nu = complex(zeros(n,1));
	invNu = complex(zeros(n,1));
	X1 = complex(zeros(n, n));
	X2 = complex(zeros(n, n));
	for i=1:n
		r = roots([1 (2 - Gamma(i)) 1]);
		Nu(i) = r(1);
		invNu(i) = r(2);
		xt1 = Zt2(:, i) - (1 / sqrt(Nu(i))) * Zt1(:, i);
		xt2 = Zt2(:, i) - sqrt(Nu(i)) * Zt1(:, i);
		X1(:, i) = xt1(1:n) + xt1(n+1:2*n);
		X2(:, i) = xt2(1:n) + xt2(n+1:2*n);
	endfor
	verbose && fprintf(': %f s\n', toc());
	
	V_II = [X1 X2];
	LAMBDA_II = [Nu; invNu];
endfunction