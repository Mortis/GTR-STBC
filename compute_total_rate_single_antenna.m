function R_sa = compute_total_rate_single_antenna(eta, sigma_w2, K, T, fixed_antenna_index)
% compute_total_rate_single_antenna: 计算固定单天线接收策略下的系统和速率

    R_sa = 0;
    for k = 1:K
        for t = 1:T
            % --- 固定单天线逻辑 ---
            % 直接获取指定天线(fixed_antenna_index)的信道增益
            % 同样，我们用eta第一行的值代表无人机天线的信道增益
            fixed_channel_gain = eta(1, fixed_antenna_index, k, t);
            
            % 计算该信道的信噪比 (SNR)
            snr = fixed_channel_gain / sigma_w2;
            
            % 累加该用户、该时隙下的速率
            R_sa = R_sa + log2(1 + snr);
        end
    end
    
    % 对时隙T进行平均
    R_sa = R_sa / T;
end