% Test of LSQR for ill-posed problems using different computing precision
% one dimensional problems

clear;
directory = pwd;
path(directory, path)
path([directory, '/regu'], path)
rng(0);  

%% test matrices
%[A, b_e, x] = shaw(1000); 
%[A,b_e,x]=deriv2(1000);  
%[A, b_e, x] = gravity(2000); 
[A, b_e, x] = heat(2000); 

[m, n] = size(A);
delta = 1e-3; 
e = randn(m, 1);  nor0 = norm(e);
e = e / norm(e) * norm(b_e) * delta;
b = b_e + e;


k = 30;
tau = 1.001; % value of \tau for discrepancy principle
[X0, R0, N0] = lsqr1(A, b, k);
[X1, R1, N1] = lsqr_s(A, b, k, "double");
[X2, R2, N2] = lsqr_s(A, b, k, "single");

er0 = zeros(k, 1); 
er1 = zeros(k, 1); 
er2 = zeros(k, 1); 

for i = 1:k
    er0(i) = norm(x - X0(:, i)) / norm(x);
    er1(i) = norm(x - X1(:, i)) / norm(x);
    er2(i) = norm(x - X2(:, i)) / norm(x);
end


[er0_opt, stop0_opt] = min(er0);
[er1_opt, stop1_opt] = min(er1);
[er2_opt, stop2_opt] = min(er1);

tol = tau * norm(e); 

for i = 1:k
    if (R0(i) < tol || R0(i) == tol)
        break;
    end
end
stop0_DP = i; er0_DP = er0(i);
for i = 1:k
    if (R0(i) < tol || R0(i) == tol)
        break;
    end
end
stop1_DP = i; er1_DP = er1(i);
for i = 1:k
    if (R0(i) < tol || R0(i) == tol)
        break;
    end
end
stop2_DP = i; er2_DP = er2(i);

[stop0_L, info0L] = corner1(R0, N0, 0); er0_L = er0(stop0_L);
[stop1_L, info1L] = corner1(R1, N1, 0); er1_L = er1(stop1_L);
[stop2_L, info2L] = corner1(R2, N2, 0); er2_L = er2(stop2_L);
[reg0_L, rho0_L,eta0_L] = l_corner(R0, N0,1:k);
[reg1_L, rho1_L,eta1_L] = l_corner(R1, N1,1:k);
[reg2_L, rho2_L,eta2_L] = l_corner(R2, N2,1:k);

%%
rel_er1 = zeros(k, 1); 
rel_er2 = zeros(k, 1); 
rel_er3 = zeros(k, 1); 

for i = 1:k
    rel_er1(i) = norm(X0(:, i) - X1(:, i)) / norm(X0(:, i));
    rel_er2(i) = norm(X0(:, i) - X2(:, i)) / norm(X0(:, i));
    rel_er3(i) = norm(X1(:, i) - X2(:, i)) / norm(X1(:, i));
end


[~, ~, B, bbeta] = LanBid(A, b, k + 1, 2);
[Q_k, R_k] = qr1(B(1:k + 1, 1:k));
d_k = diag(R_k);
D_k = diag(d_k);
Rk_til = D_k \ R_k;
cond1 = zeros(k, 1); 
cond2 = zeros(k, 1); 
cond3 = zeros(k, 1); 
bnd1 = zeros(k, 1); 
eps_s = 2^(-24);

for i = 1:k
    cond1(i) = cond(R_k(1:i, 1:i));
    cond2(i) = norm(B(1:i + 1, 1:i)) / min(d_k(1:i));
    cond3(i) = cond(Rk_til(1:i, 1:i));
    bnd1(i) = (1 + cond3(i)) * eps_s;
end

%% plot

lw = 2; l = 1:k;
figure;
semilogy(l, er0, '-bx');
hold on;
semilogy(l, er1, '-r^');
hold on;
semilogy(l, er2, '-go');
handle = legend('d', 's+d', 's+s','Location', 'northwest');
set(handle, 'Fontsize', 14);
xlabel('Iteration', 'Fontsize', 15);
ylabel('RE(k)', 'Fontsize', 15);

figure;
semilogy(rel_er1, 'b*-');
hold on;
semilogy(rel_er2, 'ro-');
handle = legend('s+d', 's+s', 'Location', 'northwest');
set(handle, 'Fontsize', 14);
xlabel('Iteration', 'Fontsize', 15);
ylabel('Relative error', 'Fontsize', 15);

figure;
semilogy(cond1, 'b*-');
hold on;
semilogy(cond2, 'rx-');
hold on;
semilogy(cond3, 'go-');
handle = legend('$\kappa(B_k)$', '$\|D_{k}^{-1}\|\cdot\|B_k\|$', '$\kappa(\widetilde{R}_k)$', 'Location', 'northeast');
set(handle, 'Fontsize', 14, 'interpreter', 'latex');
xlabel('Iteration', 'Fontsize', 15);
ylabel('Value', 'Fontsize', 15);

figure;
semilogy(rel_er3, 'b*-');
hold on;
semilogy(bnd1, 'ro-');
handle = legend('$\|\hat{x}_k-x_k\|/\|x_k\|$', 'upper bound', 'Location', 'northwest');
set(handle, 'Fontsize', 14, 'interpreter', 'latex');
xlabel('Iteration', 'Fontsize', 15);
ylabel('Relative error', 'Fontsize', 15);
