% =========================================================================
% main_vs_power_updated.m
% -------------------------------------------------------------------------
% 该脚本根据指定的物理参数（发射功率、噪声谱密度等）进行仿真。
% 它会循环遍历一个发射功率范围 (13 dBm 到 23 dBm)，
% 并最终绘制实际和速率与理论上界随功率变化的曲线。
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
% 文献值为 -40 dB，这里转换为线性尺度
tau_0_linear = 10^(-40 / 10);

% 固定系统参数
K = 8;              % 用户数
T = 10;             % 时隙数
H = 100;            % UAV 高度 (m)
alpha = 2;          % 路径损耗指数
N_antenna = 4;      % 天线数
max_iter = 10;      % 每个功率点的优化迭代次数

% 创建数组用于存储每个功率点的仿真结果
actual_rates = zeros(size(power_dbm_range));
upper_bounds = zeros(size(power_dbm_range));

% 固定用户位置和UAV初始轨迹，以确保不同功率点之间的比较是公平的
ak = rand(1,K)*500;
bk = rand(1,K)*500;
x_init = rand(1,T)*500;
y_init = rand(1,T)*500;


% --- 2. 主循环：遍历所有发射功率点 ---

fprintf('开始仿真：遍历不同发射功率...\n');
for p_idx = 1:length(power_dbm_range)
    current_power_dbm = power_dbm_range(p_idx);
    fprintf('正在运行: 发射功率 = %.1f dBm\n', current_power_dbm);

    % 将当前发射功率从dBm转换为线性尺度 (W)
    P_linear = 10^((current_power_dbm - 30) / 10);
    
    % 计算等效的 beta0，它将发射功率和参考信道增益结合在一起
    beta0 = P_linear * tau_0_linear;

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
    % 使用优化后的平均距离来计算一个平均信噪比(SNR)
    % SNR_avg = P * tau_0 * (d_avg)^-alpha / noise_power
    % SNR_avg = beta0 * (1/d_ave)^alpha / sigma_w2;
    % 根据简化的理论公式计算上界
    R_ub = K * log2(1 + SNR);
    upper_bounds(p_idx) = R_ub;
end

fprintf('所有功率点的仿真均已完成。\n');


% % --- 3. 绘图：显示最终结果 ---
% 
% figure;
% plot(power_dbm_range, actual_rates, '-o', 'LineWidth', 2, 'DisplayName', '实际和速率 (优化结果)');
% hold on;
% plot(power_dbm_range, upper_bounds, '--r', 'LineWidth', 2, 'DisplayName', '理论上界');
% grid on;
% xlabel('发射功率 (dBm)');
% ylabel('和速率 (bps/Hz)');
% title('系统和速率随发射功率的变化');
% legend('show', 'Location', 'northwest');
% --- 3. 绘图：显示最终结果 ---

% 将单位从 bps/Hz 转换为 Mbps
% 1 bps/Hz * 10e6 Hz = 10e6 bps = 10 Mbps
actual_throughput_mbps = actual_rates * (channel_bandwidth_hz / 1e6);
upper_bound_throughput_mbps = upper_bounds * (channel_bandwidth_hz / 1e6);

figure;
plot(power_dbm_range, actual_throughput_mbps, '-o', 'LineWidth', 2, 'DisplayName', '系统吞吐率');
hold on;
plot(power_dbm_range, upper_bound_throughput_mbps, '--r', 'LineWidth', 2, 'DisplayName', '理论上界');
grid on;
xlabel('发射功率 (dBm)');
ylabel('系统吞吐率 (Mbps)'); % <--- 修改Y轴标签
title('系统吞吐率随发射功率的变化'); % <--- 修改标题
legend('show', 'Location', 'northwest');