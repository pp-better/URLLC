function [ P_out ] = GetTransMatrixKrep( n, M, flag, noise, P_UE0, sinr_th, para1, para2, delta, Krep, iid, p_a_max )

P_out = zeros(n+1, n+1);
% p_a_max = 1;%M/(M+1)
for i = 0:n
    residual_n = n - i;
    if flag == 0
%         p_a = min(1, para1*M/(residual_n+1)+para2);
        p_a = min(p_a_max, para1*M/(residual_n+1)+para2);
%         if p_a >= 1
%             p_a = M/(M+1);
%         end
    else
        p_a = 1;
    end
    j_max = min(i+M, n);
           
    if p_a==1
%         p_sinr_lower = gammainc(sinr_th*(noise+P_UE0*delta^2*residual_n)/(P_UE0*(delta^2+1)),Krep*(M-residual_n+1));
        if iid == 1
            p_sinr_lower = gamcdf(sinr_th, Krep*(M-residual_n+1), (delta^2+1)*P_UE0/(noise+P_UE0*delta^2*residual_n));
        else
            p_sinr_lower = gamcdf(sinr_th, (M-residual_n+1), Krep*(delta^2+1)*P_UE0/(noise+P_UE0*delta^2*residual_n));
        end
        p_sinr_upper = 1-p_sinr_lower;
        if p_sinr_lower == 1
%             disp('L19: p_sinr_lower = 1');
            if i<n
                P_out(i+1,n+1) = 0;
            end
        else
            if p_sinr_upper == 1
                P_out(i+1, n+1) = p_sinr_upper;
                break;
            else
                if i>0 && P_out(i,i+1)==0
                    continue;
                else
                    for j = i+1:j_max
                        new_succ = j - i;
                        P_out(i+1,j+1) =  GetCombinatorialN(new_succ, residual_n)*(p_sinr_upper^new_succ)*(p_sinr_lower^(residual_n-new_succ));
                    end
                end
            end
        end
    else
        for j = i+1:j_max
            new_succ = j - i;  
            tmp = 0;
%             g_max = residual_n;
              %original one, M=4 matches well. Without
%             considerating that more than M UEs transmitting are assume to
%             fail
            g_max = min(M, residual_n);
            for g = new_succ:g_max
%                 p_sinr_lower = gammainc(sinr_th*(noise+P_UE0*delta^2*g)/(P_UE0*(delta^2+1)),Krep*(M-g+1));
                if iid == 1
                    p_sinr_lower = gamcdf(sinr_th, Krep*(M-g+1), (delta^2+1)*P_UE0/(noise+P_UE0*delta^2*g));
                else
                    p_sinr_lower = gamcdf(sinr_th, (M-g+1), Krep*(delta^2+1)*P_UE0/(noise+P_UE0*delta^2*g));
                end
                p_sinr_upper = 1-p_sinr_lower;
                if p_sinr_upper == 1
%                     disp('GetTransMatrix_g, line 46: p_sinr_upper = 1');
                    tmp = tmp + GetCombinatorialN(g, residual_n)*p_a^g*(1-p_a)^(residual_n-g);
                else
                    if p_sinr_upper > 0
                        tmp = tmp + GetCombinatorialN(g, residual_n)*GetCombinatorialN(new_succ, g)*(p_a^g)*((1-p_a)^(residual_n-g))*(p_sinr_upper^new_succ)*(p_sinr_lower^(g-new_succ));
                    end
                end
            end
            P_out(i+1,j+1) = tmp;
        end
        if P_out(i+1,j+1)<0
            disp('error: P_out<0');
        end
    end
    P_out(i+1,i+1) = 1 - sum(P_out(i+1,i+1:j_max+1));

end
if i < n
    for j = i+1:n
        P_out(j+1,n+1) = 1;
    end
end
end


