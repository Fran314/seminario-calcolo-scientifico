n = 303;
w = 1000;


A1 = complex(rand(n), rand(n))*w*w;
A0 = complex(rand(n), rand(n))*w*w;
A0 = (A0 + A0.') / 2;

[V_lqz, LAMBDA_lqz] = LQZ(A0, A1);
RR_lqz = RRes(A0, A1, V_lqz, LAMBDA_lqz);

[V_I, LAMBDA_I] = SA_I(A0, A1);
RR_I = RRes(A0, A1, V_I, LAMBDA_I);

[V_II, LAMBDA_II] = SA_II(A0, A1);
RR_II = RRes(A0, A1, V_II, LAMBDA_II);

figure(1, 'name', 'Figura 3a');
clf;
hold on;

loglog(abs(LAMBDA_lqz), RR_lqz, 'bx');
loglog(abs(LAMBDA_I), RR_I, 'ro');
loglog(abs(LAMBDA_II), RR_II, 'g^');
legend('qz', 'SA\_I', 'SA\_II');


[A0, A1] = genPQEP(3, n, w);

[V_lqz, LAMBDA_lqz] = LQZ(A0, A1);
RR_lqz = RRes(A0, A1, V_lqz, LAMBDA_lqz);

[V_I, LAMBDA_I] = SA_I(A0, A1);
RR_I = RRes(A0, A1, V_I, LAMBDA_I);

[V_II, LAMBDA_II] = SA_II(A0, A1);
RR_II = RRes(A0, A1, V_II, LAMBDA_II);

figure(2, 'name', 'Figura 3b');
clf;
hold on;

loglog(abs(LAMBDA_lqz), RR_lqz, 'bx');
loglog(abs(LAMBDA_I), RR_I, 'ro');
loglog(abs(LAMBDA_II), RR_II, 'g^');