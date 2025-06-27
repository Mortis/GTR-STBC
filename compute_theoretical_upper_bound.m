function R_upper = compute_theoretical_upper_bound(dmin, dmax, alpha, sigma_w2)
% 计算遍历速率的理论上界
%
% 输入：
%   dmin      - 最小UAV-用户距离（米）
%   dmax      - 最大UAV-用户距离（米）
%   alpha     - 路径损耗指数
%   sigma_w2  - 噪声功率
%
% 输出：
%   R_upper   - 遍历速率理论上界（bps/Hz）

    ratio = (dmax / dmin)^alpha;
    R_upper = log2(1 + 256 * ratio / sigma_w2);

    fprintf('理论上界计算：d_min = %.1f, d_max = %.1f, ratio = %.2f\n', ...
            dmin, dmax, ratio);
    fprintf('遍历速率上界 ≈ %.4f bps/Hz\n', R_upper);
end
