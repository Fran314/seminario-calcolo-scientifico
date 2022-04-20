function [K] = genKlike(m, n)
    K = zeros(n*m, n*m);
    for i=1:m
        ii = (i-1)*n + 1 : i*n;
        Kii = rand(n);
        Kii = 0.5*(Kii + Kii.');
        K(ii, ii) = Kii;
        
    end
    for i=1:m-1
        ii = (i-1)*n + 1 : i*n;
        ip1 = i*n+1 : (i+1)*n;
        Kiip1 = rand(n);
        K(ii, ip1) = Kiip1;
        K(ip1, ii) = Kiip1.';
    end
    
    Kmmp1 = rand(n);
    K(1 : n, (m-1)*n + 1 : m*n) = Kmmp1;
    K((m-1)*n + 1 : m*n, 1 : n) = Kmmp1.';
end

