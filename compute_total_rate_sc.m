function R_sc = compute_total_rate_sc(eta, sigma_w2, K, T, N_antenna)
% compute_total_rate_sc: 计算选择合并策略下的系统和速率

    R_sc = 0;
    for k = 1:K
        for t = 1:T
            % 假设用户是单天线，其信号到无人机N_antenna个接收天线的信道增益
            % 在您的eta矩阵中，任何一行的值都相同。我们用第一行代表。
            channel_gains_to_UAV_antennas = eta(1, :, k, t);
            
            % --- 选择合并核心逻辑 ---
            % 从所有接收天线中，选择信道增益最大的那一个
            best_channel_gain = max(channel_gains_to_UAV_antennas);
            
            % 计算该最佳信道的信噪比 (SNR)
            % 发射功率已经通过beta0包含在了eta中
            snr = best_channel_gain / sigma_w2;
            
            % 累加该用户、该时隙下的速率
            R_sc = R_sc + log2(1 + snr);
        end
    end
    
    % 与您原有的compute_total_rate函数一样，对时隙T进行平均
    R_sc = R_sc / T;
end