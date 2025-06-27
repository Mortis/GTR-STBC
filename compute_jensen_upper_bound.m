
function R_upper = compute_jensen_upper_bound(rho, chi, eta_avg, sigma_w2)
% 计算Jensen不等式下的遍历速率上界
% 输入：
%   rho     - [4 x K x T] 发射功率系数
%   chi     - [4 x K x T] 接收合并系数
%   eta_avg - [4 x 4 x K x T] 每条信道的平均η（通常为 β0 * d^-α）
%   sigma_w2 - 噪声功率
%
% 输出：
%   R_upper - 遍历速率上界 (bps/Hz)

[~, K, T] = size(rho);
R_upper = 0;

for k = 1:K
    for n = 1:T
        numer = 0;
        denom = 0;
        for i = 1:4
            for j = 1:4
                eta_ij = eta_avg(i,j,k,n);
                numer = numer + rho(i,k,n) * chi(j,k,n) * sqrt(sigma_w2 * eta_ij);
            end
        end
        numer = numer^2;
        
        for j = 1:4
            tmp = 0;
            for i = 1:4
                tmp = tmp + rho(i,k,n) * eta_avg(i,j,k,n);
            end
            denom = denom + sigma_w2 * chi(j,k,n)^2 * tmp;
        end
        
        R_upper = R_upper + log2(1 + numer / denom);
    end
end

R_upper = R_upper / (K * T);
end
