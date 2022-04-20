function [A0, A1] = genPQEP(m, n, omega)
    K = genKlike(m,n);
    M = genKlike(m,n);
    D = 0.2*K + 0.8*M;
    
    aA0 = complex(zeros(m*n));
    aA1 = complex(zeros(m*n));
    
    for i=1:m
        for j=1:m
            ii = (i-1)*n + 1 : i*n;
            jj = (j-1)*n + 1 : j*n;
            if(i-1 <= j && j <= i+1)
                aA0(ii, jj) = K(ii, jj) + 1i*omega*D(ii, jj) - omega*omega*M(ii,jj);
            end
            if(i == m && j == 1)
                aA1(ii, jj) = K(jj, ii) + 1i*omega*D(jj, ii) - omega*omega*M(jj, ii);
            end
        end
    end
    
    C11 = aA0(1 : n, 1 : n);
    C12 = aA0(1 : n, n+1 : (m-1)*n);
    C22 = aA0(n+1 : (m-1)*n, n+1 : (m-1)*n);
    invC22 = symmetrize(inv(C22));
    C23 = aA0(n+1 : (m-1)*n, (m-1)*n + 1 : m*n);
    C33 = aA0((m-1)*n + 1 : m*n, (m-1)*n + 1 : m*n);
    
    C = symmetrize(C11 - C12 * invC22 * C12.');
    invC = symmetrize(inv(C));
    D = - C12 * invC22 * C23;
    
    L = aA1((m-1)*n + 1 : m*n, 1 : n);
    E = symmetrize(C33 - C23.' * invC22 * C23);
    
    A1 = L * invC * D;
    A0 = symmetrize(L * invC * L.') + symmetrize(D.' * invC * D) - E;
    
    function M = symmetrize(M)
        M = 0.5 * (M + M.');
    end
end

