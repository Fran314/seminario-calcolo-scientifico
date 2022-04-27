function [REC, LAMBDA_ABS] = reciprocity(LAMBDA)
	[LAMBDA_ABS,vi]=sort(abs(LAMBDA));
	vs=LAMBDA(vi);
	REC = abs(vs.*vs(end:-1:1) - ones(size(vs)));
endfunction