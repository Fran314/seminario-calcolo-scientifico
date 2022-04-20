function [K, N, Q, Z] = rbutf(K, N)
    n = size(K, 1) / 2;

    [Q1, ~] = qr(N(1:n, 1:n));
    Q = blkdiag(Q1', eye(n));
    Z = blkdiag(eye(n), conj(Q1));
    K = Q * K * Z;
    N = Q * N * Z;
    for j=1:n-2
        for k=j+1:n-1
            % annihilate K(n+k, j) by Givens rotation in (n+k, n+k+1) plane
            G = givensl(K(n+k, j), K(n+k+1, j), 2);
            Q(n+k:n+k+1, :) = G * Q(n+k:n+k+1, :);
            K(n+k:n+k+1, :) = G * K(n+k:n+k+1, :);
            N(n+k:n+k+1, :) = G * N(n+k:n+k+1, :);
            Z(:, k:k+1) = Z(:, k:k+1) * (G.');
            K(:, k:k+1) = K(:, k:k+1) * (G.');
            N(:, k:k+1) = N(:, k:k+1) * (G.');

            % annihilate N(k+1, k) by Givens rotation in (k, k+1) plane
            G = givensl(N(k,k), N(k+1, k), 1);
            Q(k:k+1, :) = G * Q(k:k+1, :);
            K(k:k+1, :) = G * K(k:k+1, :);
            N(k:k+1, :) = G * N(k:k+1, :);
            Z(:, n+k:n+k+1) = Z(:, n+k:n+k+1) * (G.');
            K(:, n+k:n+k+1) = K(:, n+k:n+k+1) * (G.');
            N(:, n+k:n+k+1) = N(:, n+k:n+k+1) * (G.');
        endfor

        % annihilate K(2n, j) by Givens rotation in (n, 2n) plane
        G = givensl(K(n,j), K(2*n, j), 1);
        Q([n 2*n], :) = G * Q([n 2*n], :);
        K([n 2*n], :) = G * K([n 2*n], :);
        N([n 2*n], :) = G * N([n 2*n], :);
        Z(:, [n 2*n]) = Z(:, [n 2*n]) * (G');
        K(:, [n 2*n]) = K(:, [n 2*n]) * (G');
        N(:, [n 2*n]) = N(:, [n 2*n]) * (G');

        for k=n:-1:j+2
            % annihilate K(k, j) by Givens rotation in (k-1, k) plane
            G = givensl(K(k-1, j), K(k,j), 1);
            Q(k-1:k, :) = G * Q(k-1:k, :);
            K(k-1:k, :) = G * K(k-1:k, :);
            N(k-1:k, :) = G * N(k-1:k, :);
            Z(:, n+k-1:n+k) = Z(:, n+k-1:n+k) * (G.');
            K(:, n+k-1:n+k) = K(:, n+k-1:n+k) * (G.');
            N(:, n+k-1:n+k) = N(:, n+k-1:n+k) * (G.');

            % annihilate N(k, k-1) by Givens rotation in (k-1, k) plane
            G = givensr(N(k, k-1), N(k,k), 2);
            Q(n+k-1:n+k, :) = (G.') * Q(n+k-1:n+k, :);
            K(n+k-1:n+k, :) = (G.') * K(n+k-1:n+k, :);
            N(n+k-1:n+k, :) = (G.') * N(n+k-1:n+k, :);
            Z(:, k-1:k) = Z(:, k-1:k) * G;
            K(:, k-1:k) = K(:, k-1:k) * G;
            N(:, k-1:k) = N(:, k-1:k) * G;
        endfor
    endfor
    
	G = givensl(K(n,n-1), K(2*n, n-1), 1);
	Q([n 2*n], :) = G * Q([n 2*n], :);
	K([n 2*n], :) = G * K([n 2*n], :);
	N([n 2*n], :) = G * N([n 2*n], :);
	Z(:, [n 2*n]) = Z(:, [n 2*n]) * (G');
	K(:, [n 2*n]) = K(:, [n 2*n]) * (G');
	N(:, [n 2*n]) = N(:, [n 2*n]) * (G');
endfunction