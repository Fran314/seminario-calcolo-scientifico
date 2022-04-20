function G = givensr(a, b, i)
	% givensr(_, _, _) is obtained from givensl(_, _, _) 
    %   by transposing
    G = givensl(a, b, i).';
endfunction