% =========================================================================
% main_with_all_comparisons_v3.m
% -------------------------------------------------------------------------
% V3版本更新:
% 1. 在每个优化循环中加入了判断条件，确保迭代过程中的和速率是非递减的。
% 2. 将最终绘图的纵轴单位恢复为 bps/Hz (频谱效率)。
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
T = 20;             % 时隙数
H = 100;            % UAV 高度 (m)
alpha = 2;          % 路径损耗指数
N_antenna = 4;      % 天线数
max_iter = 10;       % 每个功率点的优化迭代次数

% 创建数组用于存储每个功率点的仿真结果
actual_rates = zeros(size(power_dbm_range));
upper_bounds = zeros(size(power_dbm_range));
equal_alloc_rates = zeros(size(power_dbm_range));
selection_comb_rates = zeros(size(power_dbm_range));

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
    
    % 新增：用于保证收敛非递减的变量
    rate_best_opt = -inf;
    x_best_opt = x_opt; y_best_opt = y_opt;
    rho_best_opt = rho_opt; chi_best_opt = chi_opt;

    for iter = 1:max_iter
        eta_iter = calc_eta(x_opt, y_opt, ak, bk, H, alpha, beta0, K, T);
        rho_iter = optimize_all_rho(rho_opt, chi_opt, eta_iter, sigma_w2, K, T, N_antenna);
        chi_iter = optimize_all_chi(rho_iter, chi_opt, eta_iter, sigma_w2, K, T, N_antenna);
        [x_iter, y_iter] = optimize_trajectory(x_opt, y_opt, rho_iter, chi_iter, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
        
        % 计算当前迭代的速率
        eta_current = calc_eta(x_iter, y_iter, ak, bk, H, alpha, beta0, K, T);
        rate_current = compute_total_rate(rho_iter, chi_iter, eta_current, sigma_w2, K, T, N_antenna);

        % 判断速率是否提升，确保非递减
        if rate_current >= rate_best_opt
            rate_best_opt = rate_current;
            x_best_opt = x_iter; y_best_opt = y_iter;
            rho_best_opt = rho_iter; chi_best_opt = chi_iter;
        end
        % 更新下一次迭代的输入为当前最佳状态
        x_opt = x_best_opt; y_opt = y_best_opt;
        rho_opt = rho_best_opt; chi_opt = chi_best_opt;
    end
    actual_rates(p_idx) = rate_best_opt; % 存储最终的最佳速率

    % 计算理论上界 (基于最终优化后的状态)
    [eta,d_ave] = calc_eta(x_best_opt, y_best_opt, ak, bk, H, alpha, beta0, K, T);
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
    
    rate_best_eq = -inf;
    x_best_eq = x_eq; y_best_eq = y_eq;
    
    for iter = 1:max_iter
        [x_iter, y_iter] = optimize_trajectory(x_eq, y_eq, rho_eq, chi_eq, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
        eta_current = calc_eta(x_iter, y_iter, ak, bk, H, alpha, beta0, K, T);
        rate_current = compute_total_rate(rho_eq, chi_eq, eta_current, sigma_w2, K, T, N_antenna);
        
        if rate_current >= rate_best_eq
            rate_best_eq = rate_current;
            x_best_eq = x_iter; y_best_eq = y_iter;
        end
        x_eq = x_best_eq; y_eq = y_best_eq;
    end
    equal_alloc_rates(p_idx) = rate_best_eq;

    % ==============================================================
    % C. 计算选择合并(SC)基准
    % ==============================================================
    fprintf('  正在计算选择合并基准...\n');
    objective_sc = @(qv) -compute_total_rate_sc(calc_eta(qv(1:T), qv(T+1:end), ak, bk, H, alpha, beta0, K, T), sigma_w2, K, T, N_antenna);
    q0_sc = [x_init; y_init];
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off','MaxIterations',100);
    nonlcon = @(qv) trajectory_constraints_fixed(qv, T, 50, 250, 250);
    q_opt_sc = fmincon(objective_sc, q0_sc, [], [], [], [], [], [], nonlcon, options);
    
    [eta_sc_final, ~] = calc_eta(q_opt_sc(1:T), q_opt_sc(T+1:end), ak, bk, H, alpha, beta0, K, T);
    selection_comb_rates(p_idx) = compute_total_rate_sc(eta_sc_final, sigma_w2, K, T, N_antenna);
    
    fprintf('  完成当前功率点所有计算.\n\n');
end

fprintf('所有功率点的仿真均已完成。\n');

% --- 3. 绘图：显示最终结果 (单位: bps/Hz) ---
figure;
plot(power_dbm_range, actual_rates, '-o', 'LineWidth', 2, 'DisplayName', 'GTR联合优化');
hold on;
plot(power_dbm_range, equal_alloc_rates, '-s', 'LineWidth', 2, 'DisplayName', '等功率分配');
plot(power_dbm_range, selection_comb_rates, '-^', 'LineWidth', 2, 'DisplayName', '选择合并');
plot(power_dbm_range, upper_bounds, '--r', 'LineWidth', 2, 'DisplayName', '理论上界');
grid on;
xlabel('发射功率 (dBm)');
ylabel('和速率 (bps/Hz)'); % <--- Y轴标签已改回
title('系统和速率随发射功率的变化'); % <--- 标题已修改
legend('show', 'Location', 'northwest');

fprintf('finished')