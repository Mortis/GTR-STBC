clear; clc; close all;
rng(9)
% 系统参数
K = 8;              % 用户数
T = 10;             % 时隙数
H = 100;            % UAV 高度
alpha = 2;
beta0 = 1e-4;
% sigma_w2 = 1e-5;
sigma_s = 1;
sigma_w2 = sigma_s/1e8; % 1/snr
N_antenna = 4;
d = 100;
% 有效 SNR
SNR = beta0 * (1/d)^alpha / sigma_w2;
% SNR = 1/sigma_w2;

% 换算为 dB
SNR_dB = 10 * log10(SNR);
% 参数设置
% SNR_target_dB = 0;         % 目标SNR（单位：dB）
% d_avg = 500;                % UAV平均与用户的水平距离（单位：m）
% sigma_w2 = 1e-9;            % 噪声功率
% alpha = 2.5;                % 路径损耗指数
% 
% beta0 = set_beta0_from_snr(SNR_target_dB, d_avg, alpha, sigma_w2);

% UAV 初始轨迹
x = rand(1,T)*500;
y = rand(1,T)*500;
x_=x;y_=y;
% 用户位置
ak = rand(1,K)*500;
bk = rand(1,K)*500;

% 初始化 rho 和 chi
rho = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);
chi = repmat(ones(N_antenna,1)/N_antenna, [1,K,T]);

% 主循环
max_iter = 10;
rate_history = zeros(max_iter,1);

for iter = 1:max_iter
    eta = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
    rho = optimize_all_rho(rho, chi, eta, sigma_w2, K, T, N_antenna);
    chi = optimize_all_chi(rho, chi, eta, sigma_w2, K, T, N_antenna);
    [x, y] = optimize_trajectory(x, y, rho, chi, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
    [eta,d_ave] = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
    rate_history(iter) = compute_total_rate(rho, chi, eta, sigma_w2, K, T, N_antenna);
    fprintf("Iter %d: R_avg = %.4f\n", iter, rate_history(iter));
    % figure;
    % plot(x, y, '-o', 'LineWidth', 2); hold on; scatter(ak, bk, 100, 'r', 'LineWidth', 2);plot(x_, y_, '-gx', 'LineWidth', 2);
    % xlabel('x'); ylabel('y'); title('UAV轨迹');
    % % axis equal;
    % axis([0 500 0 500]);
    % xticks(0:50:500);
    % yticks(0:50:500); 
    % legend('轨迹','用户'); grid on;
end
sigma_eta_4 = mean(eta.^2,"all");
sigma_eta_2 = mean(eta,"all");
% tol = 1e-5;              % 收敛阈值
% max_iter = 50;           % 最多迭代次数
% rate_history = zeros(max_iter, 1);
% last_rate = -inf;
% 
% for iter = 1:max_iter
%     % 重复：优化 rho, chi, trajectory
%     eta = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
%     rho = optimize_all_rho(rho, chi, eta, sigma_w2, K, T, N_antenna);
%     % [rho, rho_fval_hist] = optimize_all_rho(rho, chi, eta, sigma_w2, K, T, N_antenna);
%     % plot(rho_fval_hist{1,1}, '-o');
%     % xlabel('fmincon迭代'); ylabel('目标函数值'); title('\rho_{1,1}优化收敛');
%     chi = optimize_all_chi(rho, chi, eta, sigma_w2, K, T, N_antenna);
%     [x, y] = optimize_trajectory(rho, chi, ak, bk, H, alpha, beta0, K, T, N_antenna, sigma_w2);
%     eta = calc_eta(x, y, ak, bk, H, alpha, beta0, K, T);
%     rate = compute_total_rate(rho, chi, eta, sigma_w2, K, T, N_antenna);
% 
%     rate_history(iter) = rate;
%     fprintf("Iter %d: R_avg = %.6f\n", iter, rate);
% 
%     % 收敛判定
%     if abs(rate - last_rate) < tol
%         fprintf("已收敛（ΔR < %.1e），提前终止。\n", tol);
%         rate_history = rate_history(1:iter);
%         break;
%     end
%     last_rate = rate;
% end
SNR_ = beta0 * (1/d_ave)^alpha / sigma_w2;
SNR = sigma_eta_4/sigma_eta_2/sigma_w2;
% 绘图
figure;
plot(rate_history,'-o', 'LineWidth', 2); grid on;
xlabel('迭代'); ylabel('遍历速率 (bps/Hz)');
title('遍历速率收敛曲线');

figure;
plot(x, y, '-o', 'LineWidth', 2); hold on; scatter(ak, bk, 100, 'r', 'LineWidth', 2);plot(x_, y_, '-gx', 'LineWidth', 2);
xlabel('x'); ylabel('y'); title('UAV轨迹');
% axis equal;
axis([0 500 0 500]);
xticks(0:50:500);
yticks(0:50:500); 
legend('轨迹','用户'); grid on;

% 计算速率上界
% R_ub = compute_upper_bound(eta, sigma_w2);
dmin=100;dmax=200;
% R_ub = compute_theoretical_upper_bound(dmin, dmax, alpha, sigma_w2);
% R_ub = log2((714.14/100)^2.5*64*1+1);
R_ub = 8*log2(SNR+1);
% 绘图：遍历速率收敛曲线 + 理论上界
figure;
plot(1:max_iter, rate_history, '-o', 'LineWidth', 2); hold on;
yline(R_ub, '--r', '理论上界', 'LineWidth', 2);
xlabel('迭代次数'); ylabel('遍历速率 (bps/Hz)');
title('遍历速率收敛过程与理论上界');
legend('实际遍历速率', '理论上界');
grid on;

% rate_history=[45.5995088947400
% 47.2700824586742
% 49.0848575460691
% 50.7084960572783
% 51.9471732658736
% 52.7877127975012
% 53.2868959463574
% 53.5244980346242
% 53.6762741057644
% 53.7310462405333
% 53.7498437822193
% 53.7535670540804
% 53.7540480975915
% 53.7540989952182];
% figure;
% plot(1:14, rate_history, '-o', 'LineWidth', 2); hold on;
% xlim([1 14]);
% yline(R_ub, '--r', '理论上界', 'LineWidth', 2);
% xlabel('Iteration'); ylabel('Ergodic Rate(bps/Hz)');
% title('遍历速率收敛过程与理论上界');
% legend('Ergodic Rate', 'Upperbound');
% grid on;

% function R = compute_total_rate(rho, chi, eta, sigma_w2, K, T, N)
%     R = 0;
%     for k = 1:K
%         for t = 1:T
%             % sum_rho2_chi2 = sum(rho(:,k,t).^2) * sum(chi(:,k,t).^2);
%             % sum_chi2 = sum(chi(:,k,t).^2);
%             inner_eta = 0;
%             % for i2 = 1:N
%             %     for j2 = 1:N
%             %         inner_eta = inner_eta + rho(i2,k,t) * eta(i2,j2,k,t);
%             %     end
%             % end
%             % denom = sum_rho2_chi2 * sigma_w2 * sum_chi2 * inner_eta;
%             denom = 0;
%             for i2 = 1:N
%                 for j2 = 1:N
%                     inner_eta = inner_eta + rho(i2,k,t) * eta(i2,j2,k,t);
%                 end
%             end
%             num = 0;
%             for i = 1:N
%                 for j = 1:N
%                     num = num + rho(i,k,t) * chi(j,k,t) * sqrt(eta(i,j,k,t));
%                 end
%             end
%             num = num^2;
%             for i2 = 1:N
%                 for j2 = 1:N
%                     % denom = denom + rho(i2,k,t) * eta(i2,j2,k,t);
%                     denom = denom + rho(i,k,t).^2 * chi(j,k,t).^2 * sigma_w2 * sum(chi(:,k,t).^2) * inner_eta;
%                 end
%             end
%             % denom = sum_rho2_chi2 * sigma_w2 * sum_chi2 * inner_eta;
%             R = R + log2(1 + num / denom);
%         end
%     end
%     R = R / (K*T);
% end




