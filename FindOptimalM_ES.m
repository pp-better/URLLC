function [ M_optimal, EE_max, M_min ] = FindOptimalM_ES( N, p_arrival, K, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, pr_loss_max, delta, P_BS_tx, P_BS_C, pa_max )
% exhaustive search 

M_below_ini = max(floor(N*p_arrival),1);
M_below = M_below_ini;
% M_up = N;
M_up = 2*N;

% [ ~, p_via_mean, ~, ~, ~, ~, ~ ]= GetSystemPerformance( N, p_arrival, K, M_up, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
% while p_via_mean > pr_loss_max
%     M_up = 2*M_up;
%     [ ~, p_via_mean, ~, ~, ~, ~, ~ ]= GetSystemPerformance( N, p_arrival, K, M_up, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
% end

M_tmp_num = max(M_up-M_below_ini+1,1);

end_flag = 0;

flag = zeros(M_tmp_num, 2);

while end_flag==0
% %     p_via_old = p_via_mean;
    M = M_below + ceil((M_up-M_below)/2);
    M_id = M-M_below_ini+1;
    if flag(M_id,1) > 0
        if flag(M_id,1) > pr_loss_max
            M = M + 1;
        end
        end_flag = 1;
    else
        [ ~, p_via_mean, ~, ~, ~, ~, EE ]= GetSystemPerformance( N, p_arrival, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max);  
        flag(M_id,1) = p_via_mean;
%         if flag(M_id,1) == 0
%             flag(M_id,1) = pr_loss_max/1e5;
%         end
        if p_via_mean <= pr_loss_max
            M_up = M;
            flag(M_id,2) = EE;
        else
            M_below = M;
        end
    end
end

if end_flag == 1
    [ ~, p_via_mean, ~, ~, ~, ~, ~ ]= GetSystemPerformance( N, p_arrival, K, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max );
    if p_via_mean > pr_loss_max
        disp('FindOptimalM0:: another plan');
        M = M + 1;
    end
end
M_min = M;

M_start = M_min - M_below_ini +1;
M_ok_num = N - M_min;
for i = 0: M_ok_num
    if flag(M_start+i,2)==0
    [ ~, ~, ~, ~, ~, ~, EE ] = GetSystemPerformance( N, p_arrival, K, M_start+i+M_below_ini, pa_flag, para1, para2, payload, ber_th, tau0, noise, P_UE0, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, pa_max);
%         EE = throughput0/(P_UE_mean+(M_start+i+M_below_ini)*min(tau_mean_the/K,1)*P_BS);
        flag(M_start+i,2) = EE;
    end
end
EE_ok = flag(M_start:M_tmp_num,2);
[EE_max, M_optimal_id] = max(EE_ok); 
M_optimal = M_optimal_id + M_min - 1;

% % if M_optimal>M_max
% %     disp('FindOptimalM0:: something wrong');
% % end

end



