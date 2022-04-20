function [V_lqz, LAMBDA_lqz] = LQZ (A0, A1, verbose = false)
	verbose && fprintf('--- Linearized QZ ---\n');
	n = size(A0, 1);
	
	verbose && fprintf('- Step 1/2'); tic;
	M = [A1 zeros(n); (-A0) (-eye(n))];
	L = [zeros(n) eye(n); A1.' zeros(n)];
	verbose && fprintf(': %f s\n', toc());

	verbose && fprintf('- Step 2/2'); tic;
	[~, ~, ~, ~, V_lqz, ~, LAMBDA_lqz] = qz(M, L);
	V_lqz = V_lqz(1:n, :);
	verbose && fprintf(': %f s\n', toc());
endfunction
