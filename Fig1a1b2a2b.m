sample = 100;
%
MM_hstm = zeros(sample, 7*24);
MM_rm = zeros(sample, 7*24);
for j=1:sample
    for i=1:24
        printf('j = %d, i = %d\n', j, i);
        MM_hstm(j, (7*i - 6):(7*i)) = MaxMeanTest(50, sqrt(2)^i, true, false);
        MM_rm(j, (7*i - 6):(7*i)) = MaxMeanTest(50, sqrt(2)^i, false, false);
    endfor
endfor

figure(1, 'name', 'Figura 1a');
clf;
hold on;
%
loglog(MM_rm(1, 1:7:end), MM_rm(1, 2:7:end), 'bx-', "linewidth", 2);
loglog(MM_rm(1, 1:7:end), MM_rm(1, 3:7:end), 'bx--', "linewidth", 2);

loglog(MM_rm(1, 1:7:end), MM_rm(1, 4:7:end), 'ro-', "linewidth", 2);
loglog(MM_rm(1, 1:7:end), MM_rm(1, 5:7:end), 'ro--', "linewidth", 2);

loglog(MM_rm(1, 1:7:end), MM_rm(1, 6:7:end), 'g^-', "linewidth", 2);
loglog(MM_rm(1, 1:7:end), MM_rm(1, 7:7:end), 'g^--', "linewidth", 2);
legend('L+QZ, max', 'L+QZ, mean', 'SA\_I, max', 'SA\_I, mean', 'SA\_II, max', 'SA\_II, mean', 'location', 'northwest');
%}


figure(2, 'name', 'Figura 1b');
clf;
hold on;
%
loglog(MM_hstm(1, 1:7:end), MM_hstm(1, 2:7:end), 'bx-', "linewidth", 2);
loglog(MM_hstm(1, 1:7:end), MM_hstm(1, 3:7:end), 'bx--', "linewidth", 2);

loglog(MM_hstm(1, 1:7:end), MM_hstm(1, 4:7:end), 'ro-', "linewidth", 2);
loglog(MM_hstm(1, 1:7:end), MM_hstm(1, 5:7:end), 'ro--', "linewidth", 2);

loglog(MM_hstm(1, 1:7:end), MM_hstm(1, 6:7:end), 'g^-', "linewidth", 2);
loglog(MM_hstm(1, 1:7:end), MM_hstm(1, 7:7:end), 'g^--', "linewidth", 2);
%}

MMM_hstm = geomean(MM_hstm);
MMM_rm = geomean(MM_rm);

figure(3, 'name', 'Figura 2a');
clf;
hold on;
%
loglog(MMM_rm(1:7:end), MMM_rm(2:7:end), 'bx-', "linewidth", 2);
loglog(MMM_rm(1:7:end), MMM_rm(3:7:end), 'bx--', "linewidth", 2);

loglog(MMM_rm(1:7:end), MMM_rm(4:7:end), 'ro-', "linewidth", 2);
loglog(MMM_rm(1:7:end), MMM_rm(5:7:end), 'ro--', "linewidth", 2);

loglog(MMM_rm(1:7:end), MMM_rm(6:7:end), 'g^-', "linewidth", 2);
loglog(MMM_rm(1:7:end), MMM_rm(7:7:end), 'g^--', "linewidth", 2);
legend('L+QZ, max', 'L+QZ, mean', 'SA\_I, max', 'SA\_I, mean', 'SA\_II, max', 'SA\_II, mean', 'location', 'northwest');
%}


figure(4, 'name', 'Figura 2b');
clf;
hold on;
%
loglog(MMM_hstm(1:7:end), MMM_hstm(2:7:end), 'bx-', "linewidth", 2);
loglog(MMM_hstm(1:7:end), MMM_hstm(3:7:end), 'bx--', "linewidth", 2);

loglog(MMM_hstm(1:7:end), MMM_hstm(4:7:end), 'ro-', "linewidth", 2);
loglog(MMM_hstm(1:7:end), MMM_hstm(5:7:end), 'ro--', "linewidth", 2);

loglog(MMM_hstm(1:7:end), MMM_hstm(6:7:end), 'g^-', "linewidth", 2);
loglog(MMM_hstm(1:7:end), MMM_hstm(7:7:end), 'g^--', "linewidth", 2);

%}