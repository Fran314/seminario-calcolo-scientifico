function D = MaxMeanMinTest (n, w, hstm, verbose)
	tt = tic;
	if hstm
		[A0, A1] = genPQEP(3, n, w);
	else
		A1 = complex(rand(n), rand(n))*w*w;
		A0 = complex(rand(n), rand(n))*w*w;
		A0 = (A0 + A0.') / 2;
	endif
	verbose && fprintf('Problem generated: %f s\n\n', toc(tt)); tt = tic;

	[V_lqz, LAMBDA_lqz] = LQZ(A0, A1, verbose);
	RR_lqz = RRes(A0, A1, V_lqz, LAMBDA_lqz);
	MMM_qz = [max(RR_lqz) mean(RR_lqz) min(RR_lqz)];
	verbose && fprintf('Solved QZ: %f s\n\n', toc(tt)); tt = tic;

	[V_I, LAMBDA_I] = SA_I(A0, A1, verbose);
	RR_I = RRes(A0, A1, V_I, LAMBDA_I);
	MMM_I = [max(RR_I) mean(RR_I) min(RR_I)];
	verbose && fprintf('Solved SA_I: %f s\n\n', toc(tt)); tt = tic;

	[V_II, LAMBDA_II] = SA_II(A0, A1, verbose);
	RR_II = RRes(A0, A1, V_II, LAMBDA_II);
	MMM_II = [max(RR_II) mean(RR_II) min(RR_II)];
	verbose && fprintf('Solved SA_II: %f s\n\n', toc(tt)); tt = tic;
	
	D = [norm(A1) MMM_qz MMM_I MMM_II];
endfunction