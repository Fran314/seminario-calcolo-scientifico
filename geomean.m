function retval = geomean (input)
    sample = size(input, 1);
    retval = ones(1, size(input, 2));
    for i=1:sample
        retval .*= (input(i, :).^(1 / sample));
    endfor
endfunction