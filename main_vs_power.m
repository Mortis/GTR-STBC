% =========================================================================
% main_vs_power.m
% -------------------------------------------------------------------------
% 该脚本用于仿真遍历速率随发射功率变化的性能。
% 它会循环遍历一个发射功率范围，对每个功率点运行联合优化算法，
% 并最终绘制实际和速率与理论上界随功率变化的曲线。
% =========================================================================

clear; clc; close all;
% 设置固定的随机种子以保证结果可复现
rng(9);

% --- 1. 定义仿真参数 ---

% 定义发射功率的变化范围 (单位: dB，相对于基准功率)
power_db_range = -20:5:20;
% 将功率从dB转换为线性尺度因子
sigma_s_range = 10.^(power_db_range / 10);

% 固定系统参数
K = 8;              % 用户数
T = 10;             % 时隙数
H = 100;            % UAV 高度
alpha = 2;          % 路径损耗指数
base_beta0 = 1e-4;  % 对应基准功率(sigma_s=1)下的参考信道增益
sigma_w2 = 1e-9;    % 固定噪声功率
N_antenna = 4;      % 天线数
max_iter = 10;      % 每个功率点的优化迭代次数

% 创建数组用于存储每个功率点的仿真结果
actual_rates = zeros(size(sigma_s_range));
upper_bounds = zeros(size(sigma_s_range));

% 固定用户位置和UAV初始轨迹，以确保不同功率点之间的比较是公平的
ak = rand(1,K)*500;
bk = rand(1,K)*500;
x_init = rand(1,T)*500;
y_init = rand(1,T)*500;


% --- 2. 主循环：遍历所有发射功率点 ---

fprintf('开始仿真：遍历不同发射功率...\n');
for p_idx = 1:length(sigma_s_range)
    current_sigma_s = sigma_s_range(p_idx);
    fprintf('正在运行: 发射功率 = %.1f dB\n', power_db_range(p_idx));

    % 根据当前发射功率因子调整 beta0
    beta0 = base_beta0 * current_sigma_s;

    % 为当前功率点的优化重置所有变量
    x = x_init;
    y = y_init;
    rho = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    chi = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    rate_history_inner = zeros(max_iter,1);

    % --- 内部优化循环 (与 main_joint.m 相同) ---
    for iter = 1:max_iter
        % 计算信道增益 eta
        eta = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
        % 优化 rho
        rho = optimize_all_rho(rho, chi, eta, sigma_w2, K, T, N_antenna);
        % 优化 chi
        chi = optimize_all_chi(rho, chi, eta, sigma_w2, K, T, N_antenna);
        % 优化轨迹
        [x, y] = optimize_trajectory(x, y, rho, chi, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
        
        % 计算并记录当前迭代的和速率
        [eta,d_ave] = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
        rate_history_inner(iter) = compute_total_rate(rho, chi, eta, sigma_w2, K, T, N_antenna);
    end
    fprintf('  完成优化. 最终和速率: %.4f bps/Hz\n', rate_history_inner(end));
    sigma_eta_4 = mean(eta.^2,"all");
    sigma_eta_2 = mean(eta,"all");
    SNR_ = beta0 * (1/d_ave)^alpha / sigma_w2;
    SNR = sigma_eta_4/sigma_eta_2/sigma_w2;
    % 存储当前功率点下的最终和速率
    actual_rates(p_idx) = rate_history_inner(end);

    % 计算并存储理论上界
    % 使用优化后的平均距离来计算一个平均信噪比
    % SNR_avg = beta0 * (1/d_ave)^alpha / sigma_w2;
    % 根据原代码中的公式计算理论上界
    R_ub = K * log2(1 + SNR);
    upper_bounds(p_idx) = R_ub;
end

fprintf('所有功率点的仿真均已完成。\n');


% --- 3. 绘图：显示最终结果 ---

figure;
plot(power_db_range, actual_rates, '-o', 'LineWidth', 2, 'DisplayName', '实际和速率 (优化结果)');
hold on;
plot(power_db_range, upper_bounds, '--r', 'LineWidth', 2, 'DisplayName', '理论上界');
grid on;
xlabel('发射功率 (dB, 相对于基准)');
ylabel('遍历速率 (bps/Hz)');
title('系统和速率随发射功率的变化');
legend('show', 'Location', 'northwest');