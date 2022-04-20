function G = givensl(a, b, i)
    if(i == 1)
        G = givens(a, b);
    else
        % givensl(_, _, 2) is obtained from givensl(_, _, 1) 
        %   with a rotation of pi/2
        G = [0 (-1); 1 0] * givensl(a, b, 1);
    endif
endfunction