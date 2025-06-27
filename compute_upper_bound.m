function R_upper = compute_upper_bound(eta, sigma_w2)
% eta: [4 x 4 x K x T] 信道增益张量
% sigma_w2: 噪声功率 (scalar)

% 找到最大最小的信道增益
eta_max = max(eta(:));
eta_min = max(min(eta(:)), 1e-12);  % 避免除0

% 公式: R_upper = log2(1 + 256 * eta_max / (sigma_w2 * eta_min))
R_upper = log2(1 + (256 * eta_max) / (sigma_w2 * 1e-6));

fprintf('信道增益范围: η_min = %.2e, η_max = %.2e\n', eta_min, eta_max);
fprintf('遍历速率理论上界: R_avg <= %.4f bps/Hz\n', R_upper);
end