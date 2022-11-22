function [ path_loss ] = GetPathLoss( d_min, d_max, N, PLtype )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
d = ones(1, N);
d_mean = (d_min + d_max)/2;
if PLtype == 0
    gap = (d_max-d_min)/N;
    d_step = 0:1:N-1;
    d = d_min+gap.*d_step;
end

if PLtype == 1
    gap = (d_max-d_mean)/N;
    d(1) = d_mean;
    f = -1;
    for i = 2:N
        d(i) = d_mean + f*i*gap;
        f = f*(-1);
%         if mod(i,2)==0
%             d(i) = d_mean - floor(i/2)*gap;
%         else
%             d(i) = d_mean + floor(i/2)*gap;
%         end
    end
end

if PLtype == 2
    gap = (d_max-d_min)/N;
    d_step = 0:1:N-1;
    d = d_max-gap.*d_step;
end

path_loss = d.^(-3.76).*(10^(-3.53));%35.3+37.6lg(d)[dBm]

end

