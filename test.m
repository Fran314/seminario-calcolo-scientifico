% hstm: High Speed Train Mode (test mode)
%	true -> test with matrices in the same form as section 6
%	false -> test with random palindromic quadratic eigenvalue problem
hstm = true;
n = 100;
w = 1000;

tt = tic;
if hstm
	[A0, A1] = genPQEP(3, n, w);
else
	A1 = complex(rand(n), rand(n))*w*w;
	A0 = complex(rand(n), rand(n))*w*w;
	A0 = (A0 + A0.') / 2;
endif
fprintf('Problem generated: %f\n\n', toc(tt)); tt = tic;

clf;
hold on;

[V_lqz, LAMBDA_lqz] = LQZ(A0, A1);
RR_lqz = RRes(A0, A1, V_lqz, LAMBDA_lqz);
loglog(abs(LAMBDA_lqz), RR_lqz, 'bx');drawnow
fprintf('Solved QZ: %f\n\n', toc(tt)); tt = tic;

[V_I, LAMBDA_I] = SA_I(A0, A1);
RR_I = RRes(A0, A1, V_I, LAMBDA_I);
loglog(abs(LAMBDA_I), RR_I, 'ro');drawnow
fprintf('Solved SA_I: %f\n\n', toc(tt)); tt = tic;

[V_II, LAMBDA_II] = SA_II(A0, A1);
RR_II = RRes(A0, A1, V_II, LAMBDA_II);
loglog(abs(LAMBDA_II), RR_II, 'g^');drawnow
fprintf('Solved SA_II: %f\n\n', toc(tt)); tt = tic;

legend('qz', 'SA\_I', 'SA\_II');