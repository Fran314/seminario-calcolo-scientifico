function D = MaxMeanTest(n, w, hstm, verbose)
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
	MM_qz = [max(RR_lqz) geomean(RR_lqz)];
	verbose && fprintf('Solved QZ: %f s\n\n', toc(tt)); tt = tic;

	[V_I, LAMBDA_I] = SA_I(A0, A1, verbose);
	RR_I = RRes(A0, A1, V_I, LAMBDA_I);
	MM_I = [max(RR_I) geomean(RR_I)];
	verbose && fprintf('Solved SA_I: %f s\n\n', toc(tt)); tt = tic;

	[V_II, LAMBDA_II] = SA_II(A0, A1, verbose);
	RR_II = RRes(A0, A1, V_II, LAMBDA_II);
	MM_II = [max(RR_II) geomean(RR_II)];
	verbose && fprintf('Solved SA_II: %f s\n\n', toc(tt)); tt = tic;
	
	D = [w^2 MM_qz MM_I MM_II];
endfunction
