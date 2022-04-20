clf;

%
T = zeros(24, 10);
for i=1:24
	i
	T(i, :) = MaxMeanMinTest(50, sqrt(2)**i, false, false);
endfor
%}

T = sortrows(T);
%{
hold on;
loglog(T(:, 1), T(:, 2), 'bx-', "linewidth", 2);
loglog(T(:, 1), T(:, 3), 'bx--', "linewidth", 2);

loglog(T(:, 1), T(:, 5), 'ro-', "linewidth", 2);
loglog(T(:, 1), T(:, 6), 'ro--', "linewidth", 2);

loglog(T(:, 1), T(:, 8), 'g^-', "linewidth", 2);
loglog(T(:, 1), T(:, 9), 'g^--', "linewidth", 2);


legend('L+QZ, max', 'L+QZ, mean', 'SA\_I, max', 'SA\_I, mean', 'SA\_II, max', 'SA\_II, mean', 'location', 'northwest');
%}
