clear all;
% lambda = 0.3;%1e-3;
mu = 0.05:0.05:0.8;
mu_num = length(mu);
% N = 5;
% N0 = 15;%10:5:50;%10:10:100;10:10:150

M = 2:1:8;
M_num = length(M);


delta = 0:0.1:0.2;%
delta_num = length(delta);

K = 7;
delay_bound = 2*K;
pa_flag = 0;
para1 = 1;
para2 = 0;

tau0 = 2*1e-3/28;
payload = 20*8;
ber_th = 1e-7; 
pr_loss_max = 0.5e-5;
% P_beta = 0;
P_UE0 = 5e-13;
B = 30*30e3;
noise = 10^(-174/10)*1e-3*B;
sinr_th = 2^(sqrt(log2(exp(1))/tau0/B)*qfuncinv(ber_th)+payload/tau0/B)-1;%2^{\Omega}-1
 
unreliability_max = 1e-5;
% sum_times = 1e6;%1e7 is OK
P_beta = 0;

xi_dBm = -90;
xi = 1e-3.*10^(xi_dBm/10);


P_BS_C_dBm = 17;
P_BS_C = 10^(P_BS_C_dBm/10)*1e-3;

P_BS_tx_dBm = 30;
P_BS_tx = 10^(P_BS_tx_dBm/10)*1e-3;

pa_max = 1;

d_min = 50;
d_max = 150;
PLtype = 0;

num1 = mu_num;
num2 = M_num;
num3 = delta_num;
N_optimal = zeros(num1, num2, num3);
EE_max = zeros(num1, num2, num3);
N_max = zeros(num1, num2, num3);

for k = 1:num3
    for j = 1:num2
        for i = 1:num1
            [ N_optimal0, EE_max0, N_max0 ] = FindOptimalN_ES( M(j), mu(i), K, pa_flag, para1, para2, payload, ber_th, tau0, noise, xi, sinr_th, unreliability_max, d_max, d_min, PLtype, delta(k), P_BS_tx, P_BS_C, pa_max );
            N_optimal(i,j,k) = N_optimal0;
            EE_max(i,j,k) = EE_max0;
            N_max(i,j,k) = N_max0;
            disp(strcat('MU:',num2str(mu(i)),', M:',num2str(M(j)),', k:',num2str(k), ', EE: ', num2str(EE_max0),', N_max: ',num2str(N_max0),', N_optimal: ',num2str(N_optimal0)));
        end
    end
end
% save('data-find-N-mu-M.mat');
% disp(N_optimal);
% disp(EE_max);

[x,y]=meshgrid(M,mu);
figure (1)
for k =1:num3
    mesh(x,y,squeeze(EE_max(:,:,k)));
    hold on;
end

figure (2)
for k =1:num3
    mesh(x,y,squeeze(N_optimal(:,:,k)));
    hold on;
end

% figure (1)
% for j = 1:num2
%     plot(mu, EE_max(:,j));
%     hold on;
% end
% figure (2)
% for j = 1:num2
%     plot(mu, N_optimal(:,j));
%     hold on;
% end
% % figure (3)
% % for j = 1:num2
% %     plot(mu, N_max(:,j));
% %     hold on;
% % end