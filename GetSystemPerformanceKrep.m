function [ tau_mean, p_via_mean, throughput_mean, P_UE_mean, P_BS_C_mean, P_BS_TX_mean, EE ] = GetSystemPerformanceKrep( N, mu, delay_bound, M, pa_flag, para1, para2, payload, ber_th, tau0, noise, xi, sinr_th, path_loss, delta, P_BS_tx, P_BS_C, Krep, iid, p_a_max, Kf )
%Kmax: the maximum times of jump for the Markov chain
p_n = zeros(1, N+1);
% p_arrival = 1-exp(-lambda);
for n = 0:N
    p_n(n+1) = prod(1:N)/(prod(1:n)*prod(1:N-n))*mu^n*(1-mu)^(N-n);
end
tau_n = zeros(1, N);
p_reliable_n = zeros(1, N);
p_reliable_x = zeros(N, N-1);
P_UE_n = zeros(1, N);
P_BS_TX_n = zeros(1, N);

Kmax = floor(delay_bound/(Krep+Kf));
% p_a_max = 1;

for n = 1:N
    P_matrix = GetTransMatrixKrep( n, M, pa_flag, noise, xi, sinr_th, para1, para2, delta, Krep, iid, p_a_max );
    Q_matrix = P_matrix(1:n, 1:n);
    R_matrix = P_matrix(1:n, n+1);
%     N_matrix = inv(eye(n,n)-Q_matrix);
%     if sum(sum(isnan(N_matrix)))>0 ||sum(sum(isinf(N_matrix)))>0
%         disp('L22: something wrong!!!');
%     end
%     c = ones(n, 1);
%     z = N_matrix*c;
    
%     p_reliable_n(n) = 0;
    for k = 0:Kmax-1
        Q_tmp0 = (Q_matrix^k)*R_matrix;
        p_reliable_n(n) = p_reliable_n(n) + Q_tmp0(1);
        tau_n(n) = tau_n(n) + Q_tmp0(1)*(k+1);
    end
    Q_tmp0 = Q_matrix^Kmax;
    p_reliable_x(n,1:n-1) = Q_tmp0(1,2:n);
    
%     if z(1) > K
%         tau_n(n) = p_reliable_n(n)*z(1)+(1-p_reliable_n(n))*K;%min(z(1),K);
%     else
%         tau_n(n) = z(1);
%     end
%     tau_n(n) = min(z(1),K);
    tau_n(n) = tau_n(n)  +  (1-p_reliable_n(n))*Kmax;
    
    if pa_flag == 0
%         p_a = min(1, para1*M/(n+1)+para2);
        p_a = min(p_a_max, para1*M/(n+1)+para2);
    else
        p_a = 1;
    end
    P_UE_n(n) = n*p_a;
    P_BS_TX_n(n) = (1-p_a)^n;
    for k = 1:Kmax-1
        Q_tmp1 = Q_matrix^k;
        for j = 0:n-1
            if pa_flag == 0
%                 p_a = min(1, para1*M/(n-j+1)+para2);
                p_a = min(p_a_max, para1*M/(n-j+1)+para2);
            else
                p_a = 1;
            end
            P_UE_n(n) = P_UE_n(n) + Q_tmp1(1, j+1)*(n-j)*p_a;
            P_BS_TX_n(n) = P_BS_TX_n(n) + Q_tmp1(1, j+1)*(1-p_a)^(n-j);
        end
    end
end
tau_mean_the = sum(p_n(2:N+1).*tau_n);
tau_mean = tau_mean_the/Kmax;
% 
% for n = 1:N
% %     p_sinr_lower = gammainc(sinr_th*(noise+xi*delta^2)/(xi*(delta^2+1)),Krep*M);
%     if iid == 1
%         p_sinr_lower = gamcdf(sinr_th, Krep*M, (delta^2+1)*xi/(noise+xi*delta^2));
%     else
%         p_sinr_lower = gamcdf(sinr_th, M, Krep*(delta^2+1)*xi/(noise+xi*delta^2));
%     end
%     tmp_sinr = (1-p_sinr_lower);%^n;
%     p_reliable_n(n) = p_reliable_n(n) * tmp_sinr; 
% end
p_via_n = 1-p_reliable_n;
p_via_mean = sum(p_n(2:N+1).*p_via_n);%(1-p_n(1))*
throughput1 = sum(p_n(2:N+1).*(1:N).*payload.*p_reliable_n)/((Krep+Kf)*Kmax*tau0);
throughput0 = 0;
for n = 1:N
    tmp = 0;
    for i = 1:n-1
        tmp = tmp + p_reliable_x(n,i)*i;
    end
    throughput0 = throughput0 + p_n(n+1)*tmp;
end
throughput_mean = throughput0*payload/((Krep+Kf)*Kmax*tau0) + throughput1;
P_UE = xi./path_loss;
P_UE_mean = sum(P_UE)*sum(p_n(2:N+1).*P_UE_n)/(Kmax*N);
P_BS_C_mean = M*P_BS_C*tau_mean; 
P_BS_TX_mean = P_BS_tx*tau_mean;%*(1-sum(p_n(2:N+1).*P_BS_TX_n)/K); 
EE = throughput_mean/(P_UE_mean*Krep/(Krep+Kf)+P_BS_TX_mean/(Krep+Kf)+P_BS_C_mean);
end



