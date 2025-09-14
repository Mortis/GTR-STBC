% =========================================================================
% main_with_all_comparisons_v5.m
% -------------------------------------------------------------------------
% V5版本更新:
% 1. 为每一个发射功率点，都独立生成一张和速率随迭代次数变化的收敛曲线图。
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
max_iter = 10;      % 每个功率点的优化迭代次数

% 创建数组用于存储每个功率点的最终仿真结果
actual_rates = zeros(size(power_dbm_range));
upper_bounds = zeros(size(power_dbm_range));
equal_alloc_rates = zeros(size(power_dbm_range));
selection_comb_rates = zeros(size(power_dbm_range));

% --- 新增：使用元胞数组存储每个功率点的收敛历史 ---
convergence_hist_opt = cell(1, length(power_dbm_range));
convergence_hist_eq = cell(1, length(power_dbm_range));

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

    % 为当前收敛曲线图初始化存储向量
    rate_history_iter_opt = zeros(max_iter, 1);
    rate_history_iter_eq = zeros(max_iter, 1);

    % ==============================================================
    % A. 计算联合优化方案
    % ==============================================================
    fprintf('  正在计算联合优化方案...\n');
    x_opt = x_init; y_opt = y_init;
    rho_opt = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    chi_opt = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
    
    rate_best_opt = -inf;
    x_best_opt = x_opt; y_best_opt = y_opt;
    rho_best_opt = rho_opt; chi_best_opt = chi_opt;

    for iter = 1:max_iter
        eta_iter = calc_eta(x_opt, y_opt, ak, bk, H, alpha, beta0, K, T);
        rho_iter = optimize_all_rho(rho_opt, chi_opt, eta_iter, sigma_w2, K, T, N_antenna);
        chi_iter = optimize_all_chi(rho_iter, chi_opt, eta_iter, sigma_w2, K, T, N_antenna);
        [x_iter, y_iter] = optimize_trajectory(x_opt, y_opt, rho_iter, chi_iter, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
        
        eta_current = calc_eta(x_iter, y_iter, ak, bk, H, alpha, beta0, K, T);
        rate_current = compute_total_rate(rho_iter, chi_iter, eta_current, sigma_w2, K, T, N_antenna);

        if rate_current >= rate_best_opt
            rate_best_opt = rate_current;
            x_best_opt = x_iter; y_best_opt = y_iter;
            rho_best_opt = rho_iter; chi_best_opt = chi_iter;
        end
        x_opt = x_best_opt; y_opt = y_best_opt;
        rho_opt = rho_best_opt; chi_opt = chi_best_opt;
        
        rate_history_iter_opt(iter) = rate_best_opt; % 记录当前迭代的最佳速率
    end
    actual_rates(p_idx) = rate_best_opt;
    convergence_hist_opt{p_idx} = rate_history_iter_opt; % 存入元胞数组

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
        
        rate_history_iter_eq(iter) = rate_best_eq; % 记录当前迭代的最佳速率
    end
    equal_alloc_rates(p_idx) = rate_best_eq;
    convergence_hist_eq{p_idx} = rate_history_iter_eq; % 存入元胞数组

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

% --- 3. 绘图：显示不同方案的最终性能 (单位: bps/Hz) ---
figure;
plot(power_dbm_range, actual_rates, '-o', 'LineWidth', 2, 'DisplayName', 'GTR联合优化');
hold on;
plot(power_dbm_range, equal_alloc_rates, '-s', 'LineWidth', 2, 'DisplayName', '等功率分配');
plot(power_dbm_range, selection_comb_rates, '-^', 'LineWidth', 2, 'DisplayName', '选择合并');
plot(power_dbm_range, upper_bounds, '--r', 'LineWidth', 2, 'DisplayName', '理论上界');
grid on;
xlabel('发射功率 (dBm)');
ylabel('和速率 (bps/Hz)');
title('系统和速率随发射功率的变化');
legend('show', 'Location', 'northwest');

% --- 4. 新增：循环绘制每个功率点下的收敛曲线 ---
for p_idx = 1:length(power_dbm_range)
    figure; % 为每个功率点创建一个新的图形窗口
    plot(1:max_iter, convergence_hist_opt{p_idx}, '-o', 'LineWidth', 2, 'DisplayName', 'GTR联合优化');
    hold on;
    plot(1:max_iter, convergence_hist_eq{p_idx}, '-s', 'LineWidth', 2, 'DisplayName', '等功率分配');
    grid on;
    xlabel('迭代次数');
    ylabel('和速率 (bps/Hz)');
    % 为每个图添加动态标题
    title(sprintf('收敛曲线 (发射功率 = %d dBm)', power_dbm_range(p_idx)));
    legend('show', 'Location', 'southeast');
end