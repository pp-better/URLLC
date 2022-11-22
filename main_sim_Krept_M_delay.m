clear all;
% lambda = 0.5;%1e-3;
% p_arrival = 1-exp(-lambda);
mu = 0.5;

% N = 5;
N0 = 14;
n0_num = length(N0);

M = [2,4];%2:1:10;2, 
M_num = length(M);

Krep = 1:1:3;
Kf = 2;%:1:3;
Krep_num = length(Krep);

% K = 7;
delay_bound = 8:4:32;
delay_num = length(delay_bound);

pa_flag = 0;
para1 = 1;
para2 = 0;
% pa_num = length(para2);
% tau_mean_the = zeros(n0_num, M_num);
num1 = delay_num;
num2 = M_num;
num3 = Krep_num;
% P_UE_mean = zeros(num1, num2);
% tau_mean = zeros(num1, num2);

tau0 = 2*1e-3/28;
payload = 20*8;%20*8;
ber_th = 1e-7; 
pr_loss_max = 0.5e-5;
% P_beta = 0;
P_UE0 = 5e-13;
B = 30*30e3;%30
noise = 10^(-174/10)*1e-3*B;
sinr_th = 2^(sqrt(log2(exp(1))/tau0/B)*qfuncinv(ber_th)+payload/tau0/B)-1;%2^{\Omega}-1

% P_BS_dBm = 15:5:25;%dBm
P_BS_C_dBm = 17;
P_BS_C = 10^(P_BS_C_dBm/10)*1e-3;

P_BS_tx_dBm = 30;
P_BS_tx = 10^(P_BS_tx_dBm/10)*1e-3;

delta = 0;%0.1

pa_max = 1;

N = N0;
gap = (150-50)/N;
d_step = 0:1:N-1;
d = 50+gap.*d_step;
path_loss = d(1:N).^(-3.76).*(10^(-3.53));
P_UE = P_UE0./path_loss;


P_BS_TX_mean_theoretical = zeros(num1, num2, num3);
P_BS_C_mean_theoretical = zeros(num1, num2, num3);
P_UE_mean_theoretical = zeros(num1, num2, num3);
p_via_mean_theoretical = zeros(num1, num2, num3);
throughput_theoretical = zeros(num1, num2, num3);
EE_theoretical = zeros(num1, num2, num3);
tau_mean_theoretical = zeros(num1, num2, num3);

for i = 1:num1    
    for j = 1:num2 
        for k = 1:num3
%             if Krep(k)>1
            [ tau_mean00, p_via_mean00, throughput_mean00, P_UE_mean00, P_BS_C_mean00, P_BS_TX_mean00, EE00 ] = GetSystemPerformanceKrep( N, mu, delay_bound(i), M(j), pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, Krep(k), 1, pa_max, Kf );
%             else
%                 [ tau_mean00, p_via_mean00, throughput_mean00, P_UE_mean00, P_BS_C_mean00, P_BS_TX_mean00, EE00 ] = GetSystemPerformance( N, mu, floor(delay_bound(i)/(Krep(k)+1)), M(j), pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
%             end
            p_via_mean_theoretical(i,j,k) = p_via_mean00;
            throughput_theoretical(i,j,k) = throughput_mean00;
            P_BS_TX_mean_theoretical(i,j,k) = P_BS_TX_mean00;
            P_BS_C_mean_theoretical(i,j,k) = P_BS_C_mean00;
            P_UE_mean_theoretical(i,j,k) = P_UE_mean00;
            EE_theoretical(i,j,k) = EE00;
            tau_mean_theoretical(i,j,k) = tau_mean00;
        end
    end
end


x=delay_bound;
figure(2)
for j = 1:num2
    for k = 1:num3
        semilogy(x, squeeze(p_via_mean_theoretical(:,j,k)));
        hold on;
    end
end




