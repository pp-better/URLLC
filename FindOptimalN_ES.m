function [ N_optimal, EE_max, N_max ] = FindOptimalN_ES( M, mu, K, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, pr_loss_max, d_max, d_min, PLtype, delta, P_BS_tx, P_BS_C, pa_max )
% exhaustive search 

N_below = 1;
% M_up_ini = N;
N_up = ceil(M/mu);

% gap = (d_max-d_min)/N_up;
% d_step = 0:1:N_up-1;
% d = 50+gap.*d_step;
% path_loss = d(1:N_up).^(-3.76).*(10^(-3.53));%35.3+37.6lg(d)[dBm]
path_loss = GetPathLoss( d_min, d_max, N_up, PLtype );
[~, p_via_mean, ~, ~, ~, ~, ~ ]= GetSystemPerformance( N_up, mu, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
while p_via_mean <= pr_loss_max
    N_up = 2*N_up;
%     gap = (d_max-d_min)/N_up;
%     d_step = 0:1:N_up-1;
%     d = 50+gap.*d_step;
%     path_loss = d(1:N_up).^(-3.76).*(10^(-3.53));%35.3+37.6lg(d)[dBm]
    path_loss = GetPathLoss( d_min, d_max, N_up, PLtype );
    [ ~, p_via_mean, ~, ~, ~, ~, ~ ]= GetSystemPerformance( N_up, mu, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
end
% M_min = M_up;

N_tmp_num = max(N_up,1);

flag = zeros(N_tmp_num, 2);

% for N = 1:N_tmp_num
%     gap = (d_max-d_min)/N;
%     d_step = 0:1:N-1;
%     d = 50+gap.*d_step;
%     path_loss = d(1:N).^(-3.76).*(10^(-3.53));%35.3+37.6lg(d)[dBm]
%     [ tau_mean_the, p_via_mean, throughput0, P_UE_mean ] = GetSystemPerformance( N, p_arrival, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss);  
%     flag(N,1) = p_via_mean;
%     EE = throughput0/(P_UE_mean+M*min(tau_mean_the/K,1)*P_BS);
%     flag(M_id,2) = EE;
% end


end_flag = 0;


while end_flag==0
% %     p_via_old = p_via_mean;
    N = N_below + ceil((N_up-N_below)/2);
    
%     gap = (d_max-d_min)/N;
%     d_step = 0:1:N-1;
%     d = 50+gap.*d_step;
%     path_loss = d(1:N).^(-3.76).*(10^(-3.53));%35.3+37.6lg(d)[dBm]
    
    path_loss = GetPathLoss( d_min, d_max, N, PLtype );
    
    if flag(N,1) > 0
        if flag(N,1) > pr_loss_max
            N = N - 1;
        end
        end_flag = 1;
    else
        [~, p_via_mean, ~, ~, ~, ~, EE]= GetSystemPerformance( N, mu, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max);  
        flag(N,1) = p_via_mean;
%         if flag(M_id,1) == 0
%             flag(M_id,1) = pr_loss_max/1e5;
%         end
        if p_via_mean <= pr_loss_max
            N_below = N;
            flag(N,2) = EE;
        else
            N_up = N;
        end
    end
end

if end_flag == 1
    [ ~, p_via_mean, ~, ~, ~, ~, ~ ] = GetSystemPerformance( N, mu, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
    if p_via_mean > pr_loss_max
        disp('FindOptimalM0:: another plan');
        N = N - 1;
    end
end
N_max = N;

if N < 1
    disp('URLLC can not be supported in these system setting!');
    N_optimal = 0;
    EE_max = 0;
    N_max = 0;
    return;
end

N_end = N_max;
% M_ok_num = N - M_min;
for i = 1: N_end
    if flag(i,2)==0
%         gap = (d_max-d_min)/i;
%         d_step = 0:1:i-1;
%         d = 50+gap.*d_step;
%         path_loss = d(1:i).^(-3.76).*(10^(-3.53));%35.3+37.6lg(d)[dBm]
        path_loss = GetPathLoss( d_min, d_max, i, PLtype );

        [ ~, ~, ~, ~, ~, ~, EE ] = GetSystemPerformance( i, mu, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max);
        flag(i,2) = EE;
    end
end
EE_ok = flag(1:N_end,2);
[EE_max, N_optimal] = max(EE_ok);
% % if M_optimal>M_max
% %     disp('FindOptimalM0:: something wrong');
% % end

end



