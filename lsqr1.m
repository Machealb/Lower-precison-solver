function [X, R, N] = lsqr1(A, b, k)
    % LSQR with updating procedure, double precision

    % Haibo Li, 2020.6.10

    [~, n] = size(A);
    X = zeros(n, k); 
    R = zeros(k, 1); 
    N = zeros(k, 1); 
    [~, V, B, bbeta] = LanBid(A, b, k + 1, 2);

    % Prepare for LSQR iteration.
    w = V(:, 1);
    phi_bar = bbeta;
    rho_bar = B(1, 1);
    x = zeros(n, 1);

    for l = 1:k
        fprintf('LSQR for updating x_k: iteration %d\n', k);

        % Construct and apply orthogonal transformation.
        % rho = norm([rho_bar,B(l+1,l)]);
        rho = sqrt(rho_bar^2 + B(l + 1, l)^2);
        c1 = rho_bar / rho;
        s1 = B(l + 1, l) / rho;
        theta = s1 * B(l + 1, l + 1);
        rho_bar = -c1 * B(l + 1, l + 1);
        phi = c1 * phi_bar;
        phi_bar = s1 * phi_bar;

        % Update the solution.
        x = x + (phi / rho) * w;
        w = V(:, l + 1) - (theta / rho) * w;
        X(:, l) = x;
        R(l) = phi_bar;
        N(l) = norm(x);
    end

end
