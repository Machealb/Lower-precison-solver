
clear, clc;
directory = pwd;
path(directory, path)
path([directory, '/IRtools-master'], path)
path([directory, '/regu'], path)
path([directory, '/RestoreTools/IterativeMethods'], path)
IRtools_setup;
rng(0);  %

%% create data
% image 1
% optblur = PRset('trueImage', 'hst','BlurLevel','medium','BC', 'zero','CommitCrime', 'on');
% n = 128;
% [A, b_true, x_true, ProbInfo] = PRblurspeckle(n,optblur);

% image 2
img = imread('cameraman.tif');
x_true = imresize(img,[256,256]);
n = size(x_true,1); 
optblur.trueImage = reshape(x_true,n,n);
optblur.BlurLevel = 'severe';
optblur.CommitCrime = 'on';
optblur.BC = 'zero';
[A, b_true, x_true, ProbInfo3] = PRblurdefocus(n,optblur);

% add noise
nel = 1e-3;  
b = PRnoise(b_true, 'gauss', nel);  % Observed Image
N = n*n;

k = 100;
tau = 1.001; % value of \tau for discrepancy principle
tol = tau * norm(b-b_true);
[X0, R0, N0] = lsqr1(A, b, k);
[X1, R1, N1] = lsqr2(A, single(b), k, "double");
[X2, R2, N2] = lsqr2(A, single(b), k, "single");

er0 = zeros(k, 1); 
er1 = zeros(k, 1); 
er2 = zeros(k, 1); 
for i = 1:k
    er0(i) = norm(x_true - X0(:, i)) / norm(x_true);
    er1(i) = norm(x_true - X1(:, i)) / norm(x_true);
    er2(i) = norm(x_true - X2(:, i)) / norm(x_true);
end

[er0_opt, stop0_opt] = min(er0);  X0_opt = X0(:,stop0_opt);
[er1_opt, stop1_opt] = min(er1);  X1_opt = X1(:,stop1_opt);
[er2_opt, stop2_opt] = min(er1);  X2_opt = X0(:,stop2_opt);

for i = 1:k
    if (R0(i) <= tol )
        break;
    end
end
stop0_DP = i;  er0_DP = er0(i);  X0_DP = X0(:,stop0_DP);
for i = 1:k
    if (R0(i) <= tol )
        break;
    end
end
stop1_DP = i;  er1_DP = er1(i);  X1_DP = X1(:,stop1_DP);
for i = 1:k
    if (R0(i) <= tol)
        break;
    end
end
stop2_DP = i;  er2_DP = er2(i);  X2_DP = X2(:,stop2_DP);

[stop0_L, info0L] = corner1(R0, N0, 0);  er0_L = er0(stop0_L);
[stop1_L, info1L] = corner1(R1, N1, 0);  er1_L = er1(stop1_L);
[stop2_L, info2L] = corner1(R2, N2, 0);  er2_L = er2(stop2_L);
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

%%
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

%%
lw = 2; l = 1:k;
figure;
semilogy(l, er0, '-bx');
hold on;
semilogy(l, er1, '-r^');
hold on;
semilogy(l, er2, '-go');
handle = legend('d', 's+d', 's+s', 'Location', 'northeast');
set(handle, 'Fontsize', 14);
xlabel('Iteration', 'Fontsize', 15);
ylabel('RE(k)', 'Fontsize', 15);

figure;
semilogy(rel_er1, 'b*-');
hold on;
semilogy(rel_er2, 'rx-');
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
handle = legend('$\kappa(B_k)$', '$\|D_{k}^{-1}\|\cdot\|B_k\|$', '$\kappa(\widetilde{R}_k)$', 'Location', 'northwest');
set(handle, 'Fontsize', 14,'interpreter', 'latex');
xlabel('Iteration', 'Fontsize', 15);
ylabel('Value', 'Fontsize', 15);

figure;
semilogy(rel_er3, 'b*-');
hold on;
semilogy(bnd1, 'rx-');
handle = legend('$\|\hat{x}_k-x_k\|/\|x_k\|$', 'upper bound', 'Location', 'northwest');
set(handle, 'Fontsize', 14, 'interpreter', 'latex');
xlabel('Iteration', 'Fontsize', 15);
ylabel('Relative error', 'Fontsize', 15);


figure;
subplot(1,1,1);
imshow(reshape(x_true(:), n, n), []); title('true image', 'Fontsize', 16);

figure;
subplot(1,1,1);
img1 = imresize(reshape(b(:), n, n),[360,360]);
imshow(img1,[]); title('blurred image', 'Fontsize', 16);

figure;
subplot(1,1,1);
img2 = imresize(reshape(X0_opt(:), n, n),[360,360]);
imshow(img2,[]); title(sprintf('k=%d', stop0_opt), 'Fontsize', 16);


% semilogy(l,err0, '-*k','LineWidth',lw,'MarkerIndices',1:10:k,'MarkerSize',4), hold on;
% semilogy(l,err1, '-.b','LineWidth',lw)
% semilogy(l,err2, '--r','LineWidth',lw)
% semilogy(l,err3, '-c','LineWidth',lw)
% semilogy(stop0_opt, er0_opt,'bV', 'MarkerSize',16, 'LineWidth',lw)
% semilogy(stop1_opt, er1_opt,'mx', 'MarkerSize',16, 'LineWidth',lw)
% semilogy(stop2, err2(stop2),'gD', 'MarkerSize',16, 'LineWidth',lw)
% semilogy(stop3, err3(stop3),'ro', 'MarkerSize',16, 'LineWidth',lw)
% xlabel('Iteration'), ylabel('$RE(k)$','Interpreter','latex')
% legend('LSQR', 'JBDQR','JBD-WGCV','JBD-sec') 
% axis([0,k,0,1])


