function [X, R, N] = lsqr2(A, b, k, precision)
    % LSQR with updating procedure (another version):
    % the Lanczos bidiagonalization is implemented using single precision;
    % the updating procedure's computing precision can be chosen

    % Haibo Li, 2022.6.12

    [~, n] = size(A);
    X = zeros(n, k);
    R = zeros(k, 1); 
    N = zeros(k, 1); 
    [~, V, B, bbeta] = LanBid(A, single(b), k + 1, 2);
    B = double(B); bbeta = double(bbeta);

    if (precision == "double")
        V = double(V);
    end

    % Prepare for LSQR iteration.
    w = V(:, 1);
    phi_bar = bbeta;
    rho_bar = B(1, 1);
    x = zeros(n, 1);

    if (precision == "single")
        x = single(x);
    end
    
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
        if (precision == "double")
            x = x + (phi / rho) * w;
            w = V(:, l + 1) - (theta / rho) * w;
            X(:, l) = x;
            R(l) = phi_bar;
            N(l) = norm(x);
        elseif (precision == "single")
            x = x + single(phi / rho) * w;
            w = V(:, l + 1) - single(theta / rho) * w;
            X(:, l) = x;
            R(l) = phi_bar;
            N(l) = norm(x);
        end

    end

    X = double(X); R = double(R); N = double(N);
end
