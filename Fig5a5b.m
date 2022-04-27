n = 303;
w = 1000;


[A0, A1] = genPQEP(3, n, w);

[~, LAMBDA_lqz] = LQZ(A0, A1);
[REC_lqz, LAMBDA_ABS_lqz] = reciprocity(LAMBDA_lqz);

[~, LAMBDA_I] = SA_I(A0, A1);
[REC_I, LAMBDA_ABS_I] = reciprocity(LAMBDA_I);

[~, LAMBDA_II] = SA_II(A0, A1);
[REC_II, LAMBDA_ABS_II] = reciprocity(LAMBDA_II);


figure(1, 'name', 'Figura 5a');
clf;
hold on;

loglog(LAMBDA_ABS_lqz, REC_lqz, 'bx');
legend('qz');


figure(2, 'name', 'Figura 5b');
clf;
hold on;

loglog(LAMBDA_ABS_I, REC_I, 'ro');
loglog(LAMBDA_ABS_II, REC_II, 'g^');
legend('SA\_I', 'SA\_II');