% =========================================================================
% main_with_all_comparisons.m
% -------------------------------------------------------------------------
% 该脚本在之前版本的基础上，新增了“选择合并”作为第二个基准对比。
% =========================================================================
clear; clc; close all;
% 设置固定的随机种子以保证结果可复现
rng(9);

% --- 1. 定义仿真参数 (根据文献设定) ---
% 定义发射功率的变化范围 (单位: dBm)
power_dbm_range = -24:5:23;

% 噪声参数
noise_psd_dbm_hz = -160; % 噪声功率谱密度 (dBm/Hz)
channel_bandwidth_hz = 10e6; % 信道带宽 (10 MHz)
% 计算总噪声功率 (线性尺度)
noise_power_dbm = noise_psd_dbm_hz + 10 * log10(channel_bandwidth_hz);
sigma_w2 = 10^((noise_power_dbm - 30) / 10); % 转换为瓦特(W)

% 信道参考功率 (对应文献中的 tau_0)
tau_0_linear = 10^(-40 / 10);

% 固定系统参数
K = 8;              % 用户数
T = 10;             % 时隙数
H = 100;            % UAV 高度 (m)
alpha = 2;          % 路径损耗指数
N_antenna = 4;      % 天线数
max_iter = 3;      % 每个功率点的优化迭代次数

% 创建数组用于存储每个功率点的仿真结果
actual_rates = zeros(size(power_dbm_range));
upper_bounds = zeros(size(power_dbm_range));
equal_alloc_rates = zeros(size(power_dbm_range));
selection_comb_rates = zeros(size(power_dbm_range)); % <-- 新增：存储选择合并结果

% 固定用户位置和UAV初始轨迹
ak = rand(1,K)*500;
bk = rand(1,K)*500;
x_init = rand(1,T)*500;
y_init = rand(1,T)*500;


% --- 2. 主循环：遍历所有发射功率点 ---
fprintf('开始仿真：遍历不同发射功率...\n');
for p_idx = 1:length(power_dbm_range)
    current_power_dbm = power_dbm_range(p_idx);
    fprintf('正在运行: 发射功率 = %.1f dBm\n', current_power_dbm);
    
    P_linear = 10^((current_power_dbm - 30) / 10);
    beta0 = P_linear * tau_0_linear;

    % ==============================================================
    % A. 计算联合优化方案
    % ==============================================================
    fprintf('  正在计算联合优化方案...\n');
    x_opt = x_init; y_opt = y_init;
    rho_opt = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    chi_opt = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    for iter = 1:max_iter
        eta = calc_eta(x_opt, y_opt, ak, bk, H, alpha, beta0, K, T);
        rho_opt = optimize_all_rho(rho_opt, chi_opt, eta, sigma_w2, K, T, N_antenna);
        chi_opt = optimize_all_chi(rho_opt, chi_opt, eta, sigma_w2, K, T, N_antenna);
        [x_opt, y_opt] = optimize_trajectory(x_opt, y_opt, rho_opt, chi_opt, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
    end
    [eta,d_ave] = calc_eta(x_opt, y_opt, ak, bk, H, alpha, beta0, K, T);
    actual_rates(p_idx) = compute_total_rate(rho_opt, chi_opt, eta, sigma_w2, K, T, N_antenna);
    sigma_eta_4 = mean(eta.^2,"all"); sigma_eta_2 = mean(eta,"all");
    SNR = sigma_eta_4/sigma_eta_2/sigma_w2;
    upper_bounds(p_idx) = K * log2(1 + SNR);
    
    % ==============================================================
    % B. 计算等功率分配基准
    % ==============================================================
    fprintf('  正在计算等功率分配基准...\n');
    x_eq = x_init; y_eq = y_init;
    rho_eq = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    chi_eq = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    for iter = 1:max_iter
        [x_eq, y_eq] = optimize_trajectory(x_eq, y_eq, rho_eq, chi_eq, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
    end
    [eta_eq, ~] = calc_eta(x_eq, y_eq, ak, bk, H, alpha, beta0, K, T);
    equal_alloc_rates(p_idx) = compute_total_rate(rho_eq, chi_eq, eta_eq, sigma_w2, K, T, N_antenna);

    % ==============================================================
    % C. 新增：计算选择合并(SC)基准
    % ==============================================================
    fprintf('  正在计算选择合并基准...\n');
    x_sc = x_init; y_sc = y_init;
    
    % --- 为选择合并策略独立优化轨迹 ---
    % 定义目标函数：最大化SC速率，即最小化其负数
    objective_sc = @(qv) -compute_total_rate_sc(calc_eta(qv(1:T), qv(T+1:end), ak, bk, H, alpha, beta0, K, T), sigma_w2, K, T, N_antenna);
    
    % 使用与 optimize_trajectory.m 内部相同的 fmincon 设置
    q0_sc = [x_sc; y_sc];
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','MaxIterations',100);
    nonlcon = @(qv) trajectory_constraints_fixed(qv, T, 50, 250, 250); % 使用相同的约束

    % 直接调用fmincon进行轨迹优化
    q_opt_sc = fmincon(objective_sc, q0_sc, [], [], [], [], [], [], nonlcon, options);
    
    x_sc_opt = q_opt_sc(1:T);
    y_sc_opt = q_opt_sc(T+1:end);
    
    % 使用优化后的轨迹计算最终的SC速率
    [eta_sc_final, ~] = calc_eta(x_sc_opt, y_sc_opt, ak, bk, H, alpha, beta0, K, T);
    selection_comb_rates(p_idx) = compute_total_rate_sc(eta_sc_final, sigma_w2, K, T, N_antenna);

    fprintf('  完成当前功率点所有计算.\n\n');
end

fprintf('所有功率点的仿真均已完成。\n');

% --- 3. 绘图：显示最终结果 ---
% 将单位从 bps/Hz 转换为 Mbps
actual_throughput_mbps = actual_rates * (channel_bandwidth_hz / 1e6);
upper_bound_throughput_mbps = upper_bounds * (channel_bandwidth_hz / 1e6);
equal_alloc_throughput_mbps = equal_alloc_rates * (channel_bandwidth_hz / 1e6);
selection_comb_throughput_mbps = selection_comb_rates * (channel_bandwidth_hz / 1e6); %<-- 新增

figure;
plot(power_dbm_range, actual_throughput_mbps, '-o', 'LineWidth', 2, 'DisplayName', '联合优化');
hold on;
plot(power_dbm_range, equal_alloc_throughput_mbps, '-s', 'LineWidth', 2, 'DisplayName', '等功率分配');
plot(power_dbm_range, selection_comb_throughput_mbps, '-^', 'LineWidth', 2, 'DisplayName', '选择合并'); %<-- 新增对比曲线
plot(power_dbm_range, upper_bound_throughput_mbps, '--r', 'LineWidth', 2, 'DisplayName', '理论上界');
grid on;
xlabel('发射功率 (dBm)');
ylabel('系统吞吐率 (Mbps)');
title('系统吞吐率随发射功率的变化');
legend('show', 'Location', 'northwest');

