clear; clc; close all;
% rng(13)
% 系统参数
K = 8;              % 用户数
T = 10;             % 时隙数
H = 100;            % UAV 高度
alpha = 2;
beta0 = 1e-4;
sigma_s = 1;
sigma_w2 = sigma_s/1e8; % 1/snr
N_antenna = 4;

% UAV 初始轨迹
x = rand(1,T)*500;
y = rand(1,T)*500;
x_ = x; y_ = y;

% 用户位置
ak = rand(1,K)*500;
bk = rand(1,K)*500;

% 初始化 rho 和 chi
rho = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
chi = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);

% 主循环
max_iter = 6;
rate_history = zeros(max_iter,1);

% =========================================================================
% --- 新增：用于保证收敛非递减的变量 ---
rate_best = -inf;       % 用于存储迭代过程中的最佳速率
x_best = x;             % 存储最佳速率对应的x轨迹
y_best = y;             % 存储最佳速率对应的y轨迹
rho_best = rho;         % 存储最佳速率对应的rho
chi_best = chi;         % 存储最佳速率对应的chi
% =========================================================================

for iter = 1:max_iter
    % 使用上一轮的最佳状态作为本次迭代的输入
    eta_iter = calc_eta(x_best, y_best, ak, bk, H, alpha, beta0, K, T);
    rho_iter = optimize_all_rho(rho_best, chi_best, eta_iter, sigma_w2, K, T, N_antenna);
    chi_iter = optimize_all_chi(rho_iter, chi_best, eta_iter, sigma_w2, K, T, N_antenna);
    [x_iter, y_iter] = optimize_trajectory(x_best, y_best, rho_iter, chi_iter, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
    
    % 计算当前迭代产生的速率
    [eta_current, ~] = calc_eta(x_iter, y_iter, ak, bk, H, alpha, beta0, K, T);
    rate_current = compute_total_rate(rho_iter, chi_iter, eta_current, sigma_w2, K, T, N_antenna);
    
    % --- 判断速率是否提升，确保非递减 ---
    if rate_current >= rate_best
        % 如果速率有提升或保持不变，则更新最佳状态
        rate_best = rate_current;
        x_best = x_iter;
        y_best = y_iter;
        rho_best = rho_iter;
        chi_best = chi_iter;
    end
    
    % 记录每次迭代后的最佳速率，以保证rate_history是非递减的
    rate_history(iter) = rate_best;
    
    fprintf("Iter %d: R_avg = %.4f (Current Best = %.4f)\n", iter, rate_current, rate_best);
end

% 循环结束后，最终的轨迹和SNR等计算应基于最佳状态
[eta_final, d_ave] = calc_eta(x_best, y_best, ak, bk, H, alpha, beta0, K, T);
sigma_eta_4 = mean(eta_final.^2,"all");
sigma_eta_2 = mean(eta_final,"all");

SNR_ = beta0 * (1/d_ave)^alpha / sigma_w2;
SNR = sigma_eta_4/sigma_eta_2/sigma_w2;

% 绘图
figure;
plot(rate_history,'-o', 'LineWidth', 2); grid on;
xlabel('迭代'); ylabel('遍历速率 (bps/Hz)');
title('遍历速率收敛曲线 (非递减)');

figure;
plot(x_best, y_best, '-o', 'LineWidth', 2); hold on; % 使用x_best, y_best绘图
scatter(ak, bk, 100, 'r', 'LineWidth', 2);
plot(x_, y_, '-gx', 'LineWidth', 2);
xlabel('x'); ylabel('y'); title('UAV轨迹');
axis([0 500 0 500]);
xticks(0:50:500);
yticks(0:50:500); 
legend('优化后轨迹','用户','初始轨迹'); grid on;

% 计算速率上界
R_ub = 8*log2(SNR+1);

% 绘图：遍历速率收敛曲线 + 理论上界
figure;
plot(1:max_iter, rate_history, '-o', 'LineWidth', 2); hold on;
yline(R_ub, '--r', '理论上界', 'LineWidth', 2);
xlabel('迭代次数'); ylabel('遍历速率 (bps/Hz)');
title('遍历速率收敛过程与理论上界');
legend('实际遍历速率', '理论上界');
grid on;